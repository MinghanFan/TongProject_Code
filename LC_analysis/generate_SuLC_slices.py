#!/usr/bin/env python3
"""Find every slice containing Su et al.'s LC core mesh and save outline-only overlays.

Designed for the TongProject / LeonardoLupori brainAlignment directory structure.
This is a Python implementation of the established Su mesh -> QuickNII plane
workflow: it reads the mouse's original VisuAlign JSON, intersects the
actual 3-D Su core mesh with each QuickNII plane, and applies VisuAlign's
piecewise-linear marker mapping to place LC outlines on histology thumbnails.
No CCFv3 region labels or boundaries are drawn, and MATLAB is not needed.

Dependencies: numpy, scipy, pillow
    python3 -m pip install numpy scipy pillow

Change MOUSE_ID below, then run:
    python3 generate_SuLC_slices.py

Outputs (only LC-positive images):
    LC_analysis/CC3A/SuLC_only/CC3A_007_A3_1_SuLCOverlay.png
    LC_analysis/CC3A/SuLC_only/CC3A_SuLC_all_slices.csv
    LC_analysis/CC3A/SuLC_only/CC3A_SuLC_positive_slices.csv

IMPORTANT: Check a known, previously verified slice against your MATLAB
output. SciPy's triangulation and direct full-resolution piecewise-linear
interpolation need not be pixel-identical to MATLAB's 0.2x displacement-field
sampling/resizing, despite using the same control markers and mapping direction.
"""

from __future__ import annotations

import csv
import json
import math
import re
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.spatial import Delaunay, QhullError


# ================================ SETTINGS =================================
MOUSE_ID = "CC4B"  # Change only this to process a different mouse.
PROJECT_FOLDER = Path.home() / "Desktop" / "TongProject" / "TongProject_Code"

MESH_RELATIVE_PATH = Path("atlas") / "Su_LC" / "core_mesh.obj"
VOXEL_SIZE_UM = 25.0
MIDLINE_X_UM = 5700.0
# Original annotationVolume.mat dimensions (BEFORE MATLAB permute [3 1 2]).
# Only dimensions are used for Su->QuickNII coordinate conversion; the atlas
# labels and the MATLAB .mat file itself are NOT read or displayed.
ORIGINAL_ATLAS_SHAPE = (528, 320, 456)

MIN_SU_PIXELS = 1  # One or more rasterized pixels in the unwarped LC section.
SKIP_INVALID_SLICES = True
OUTLINE_COLOR = (0, 255, 255)  # Cyan; matches the prior Su LC overlays.
OUTLINE_THICKNESS = 2         # Pixels; set to 1 for a thin outline.
# ============================================================================


def _require(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{description} not found: {path}")


def _local_tag(name: str) -> str:
    return name.rsplit("}", 1)[-1]


def load_sections(info_xml: Path) -> list[dict]:
    """Read MATLAB writestruct-produced *-info.xml without guessing filenames."""
    root = ET.parse(info_xml).getroot()
    result = []
    for el in root.iter():
        if _local_tag(el.tag).lower() != "slices":
            continue
        data = {_local_tag(c.tag).lower(): (c.text or "").strip()
                for c in el if len(c) == 0}
        if not data.get("name"):
            continue  # enclosing <slices> wrapper, not an actual entry
        v = data.get("valid", "1").strip().lower()
        valid = v not in {"0", "false", "no"}
        result.append({"name": data["name"],
                       "number": data.get("number", ""), "valid": valid})
    if not result:
        raise ValueError(f"No <slices><name> entries found in {info_xml}")
    return result


def load_alignment_records(alignment_json: Path) -> dict[str, dict]:
    """Index VisuAlign records by full slide ID, not loose substring matching."""
    # JSON exported by some versions of VisuAlign contains Windows paths.
    raw = alignment_json.read_text(encoding="utf-8-sig")
    try:
        doc = json.loads(raw)
    except json.JSONDecodeError:
        # MATLAB's existing helper normalizes path separators before decode.
        doc = json.loads(raw.replace("\\", "/"))
    records = doc.get("slices", [])
    if not isinstance(records, list) or not records:
        raise ValueError(f"Missing nonempty 'slices' array: {alignment_json}")
    result = {}
    for rec in records:
        # A filename can look like '/path/CC3A_007_A3_1-thumb.png'.
        stem = Path(str(rec.get("filename", "")).replace("\\", "/")).stem
        name = re.sub(r"-thumb$", "", stem, flags=re.IGNORECASE)
        if name in result:
            raise ValueError(f"Duplicate original alignment for {name}")
        result[name] = rec
    return result


def load_triangular_obj(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read vertices and triangular faces (including OBJ 'v/vt/vn' syntax)."""
    vertices, faces = [], []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("v "):
                coords = line.split()
                vertices.append([float(x) for x in coords[1:4]])
            elif line.startswith("f "):
                parts = line.split()[1:]
                if len(parts) != 3:
                    raise ValueError("Expected triangulated Su core_mesh.obj (3 vertices/face)")
                faces.append([int(item.split("/")[0]) - 1 for item in parts])
    v = np.asarray(vertices, dtype=np.float64).reshape(-1, 3)
    f = np.asarray(faces, dtype=np.int64).reshape(-1, 3)
    if not len(v) or not len(f) or not np.isfinite(v).all():
        raise ValueError(f"Mesh contains no valid triangles: {path}")
    if f.min() < 0 or f.max() >= len(v):
        raise ValueError("Invalid OBJ triangle indices")
    return v, f


def load_bilateral_quicknii_mesh(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Replicate the verified Su xyz -> project QuickNII conversion.

    Matches su_lc_tools.m: mirrored left/right around x=5700 um and
    A=flip(flip(permute(annotationVolume,[3,1,2]),2),3).
    """
    v, faces = load_triangular_obj(path)
    mirrored = v.copy()
    mirrored[:, 0] = 2 * MIDLINE_X_UM - mirrored[:, 0]
    all_vertices = np.concatenate([v, mirrored])
    all_faces = np.concatenate([faces, faces + len(v)])

    # After permute([3,1,2]), the atlas has dimensions (456, 528, 320).
    atlas_shape = (ORIGINAL_ATLAS_SHAPE[2],
                   ORIGINAL_ATLAS_SHAPE[0], ORIGINAL_ATLAS_SHAPE[1])
    q = np.column_stack((
        1.0 + all_vertices[:, 0] / VOXEL_SIZE_UM,
        atlas_shape[1] - 1.0 - all_vertices[:, 1] / VOXEL_SIZE_UM,
        atlas_shape[2] - 1.0 - all_vertices[:, 2] / VOXEL_SIZE_UM,
    ))
    return q, all_faces


def intersect_mesh_with_plane(q: np.ndarray, faces: np.ndarray,
                              anchoring: np.ndarray, w: int, h: int) -> np.ndarray:
    """Return 2-D surface segments [x0,y0,x1,y1] in 0-based thumbnail pixels.

    Exact triangle/plane section; unlike projecting entire 3-D mesh onto 2-D.
    Formula matches the earlier MATLAB su_lc_tools('sliceMask',...).
    """
    o, u, v = anchoring[:3], anchoring[3:6], anchoring[6:9]
    n = np.cross(u, v)
    if np.linalg.norm(n) < 1e-12:
        raise ValueError("QuickNII o/u/v has degenerate plane vectors")
    n = n / np.linalg.norm(n)
    distances = (q - o) @ n
    td = distances[faces]
    crossing = (td.min(axis=1) <= 0) & (td.max(axis=1) > 0)
    fc = faces[crossing]
    if not len(fc):
        return np.empty((0, 4), dtype=np.float64)

    # Least-squares projection into the (u,v) basis of the oblique plane.
    basis = np.column_stack((u, v))  # 3x2
    backproject = np.linalg.pinv(basis)  # 2x3
    entries = []  # (triangle index, x, y)
    for edge in range(3):
        a = fc[:, edge]
        b = fc[:, (edge + 1) % 3]
        da, db = distances[a], distances[b]
        cross = ((da <= 0) & (db > 0)) | ((db <= 0) & (da > 0))
        if not cross.any():
            continue
        t = da[cross] / (da[cross] - db[cross])
        points = q[a[cross]] + t[:, None] * (q[b[cross]] - q[a[cross]])
        xy = (points - o) @ backproject.T
        # MATLAB pixel coords are [xy*W+1, xy*H+1].  Subtract 1 for Python.
        xs, ys = xy[:, 0] * w, xy[:, 1] * h
        for idx, x, y in zip(np.nonzero(cross)[0], xs, ys):
            entries.append((int(idx), float(x), float(y)))
    if not entries:
        return np.empty((0, 4), dtype=np.float64)

    # Robustly pair exactly two crossings from the same triangle.
    per_triangle: dict[int, list[tuple[float, float]]] = {}
    for i, x, y in entries:
        per_triangle.setdefault(i, []).append((x, y))
    result = []
    for pair in per_triangle.values():
        if len(pair) == 2:
            (x0, y0), (x1, y1) = pair
            if math.hypot(x1 - x0, y1 - y0) > 1e-8:
                result.append((x0, y0, x1, y1))
    return np.array(result, dtype=np.float64).reshape(-1, 4)


def rasterize_section(segments: np.ndarray, h: int, w: int) -> np.ndarray:
    """Scanline even-odd fill, following the prior su_lc_tools MATLAB code."""
    mask = np.zeros((h, w), dtype=bool)
    if not len(segments):
        return mask
    lo = max(0, math.ceil(min(segments[:, 1].min(), segments[:, 3].min())))
    hi = min(h - 1, math.floor(max(segments[:, 1].max(), segments[:, 3].max())))
    for row in range(lo, hi + 1):
        y0, y1 = segments[:, 1], segments[:, 3]
        cross = ((y0 <= row) & (y1 > row)) | ((y1 <= row) & (y0 > row))
        if not cross.any():
            continue
        s = segments[cross]
        xs = s[:, 0] + (row - s[:, 1]) * (s[:, 2] - s[:, 0]) / (s[:, 3] - s[:, 1])
        xs.sort()
        if len(xs) % 2:
            warnings.warn(f"Odd mesh intersections in scanline {row}; skipping it")
            continue
        for a, b in xs.reshape(-1, 2):
            first, last = max(0, math.ceil(a)), min(w - 1, math.floor(b))
            if first <= last:
                mask[row, first:last + 1] = True
    return mask


def warp_mask_with_original_markers(source_mask: np.ndarray,
                                  markers: np.ndarray) -> np.ndarray:
    """Apply original VisuAlign's in-plane piecewise-affine mapping (inverse warp).

    Markers: [atlas/reference x,y, histology x,y]. At each output histology
    pixel, interpolate markers' source/reference x,y and sample source mask.
    Include the four zero-displacement corners used by triangulateSlice.m.
    Computes only a safely expanded LC region for speed and uses nearest-neighbor
    sampling, as the earlier MATLAB image masks did.
    """
    h, w = source_mask.shape
    if markers.ndim != 2 or markers.shape[1] != 4 or not len(markers):
        raise ValueError("Original VisuAlign alignment requires nonempty Nx4 markers")
    if not np.isfinite(markers).all():
        raise ValueError("VisuAlign markers contain nonfinite coordinates")
    xy = np.array([
        [-0.1 * w, -0.1 * h], [1.1 * w, -0.1 * h],
        [-0.1 * w, 1.1 * h], [1.1 * w, 1.1 * h],
    ], dtype=float)
    xy = np.concatenate((xy, markers[:, 2:4]), axis=0)
    src = np.concatenate((xy[:4], markers[:, :2]), axis=0)
    try:
        triangulation = Delaunay(xy)
    except QhullError as exc:
        raise ValueError(f"Cannot triangulate VisuAlign markers: {exc}") from exc

    rows, cols = np.nonzero(source_mask)
    if not len(rows):
        return np.zeros_like(source_mask)
    # Any barycentric interpolation lies between its vertices' displacements.
    delta = src - xy
    pad_x = int(np.ceil(np.max(np.abs(delta[:, 0])))) + 4
    pad_y = int(np.ceil(np.max(np.abs(delta[:, 1])))) + 4
    x_start, x_end = max(0, int(cols.min()) - pad_x), min(w, int(cols.max()) + pad_x + 1)
    y_start, y_end = max(0, int(rows.min()) - pad_y), min(h, int(rows.max()) + pad_y + 1)

    warped = np.zeros_like(source_mask)
    for top in range(y_start, y_end, 128):
        bottom = min(y_end, top + 128)
        yy, xx = np.mgrid[top:bottom, x_start:x_end]
        # MATLAB's getDisplacement receives pixel coordinates starting at 1.
        positions = np.column_stack((xx.ravel() + 1., yy.ravel() + 1.))
        simplex = triangulation.find_simplex(positions)
        sampled_xy = positions.copy()  # Identity for rare queries outside hull.
        inside = simplex >= 0
        if inside.any():
            ids = simplex[inside]
            transform = triangulation.transform[ids]
            b = np.einsum('nij,nj->ni', transform[:, :2],
                          positions[inside] - transform[:, 2])
            weights = np.column_stack((b, 1 - b.sum(axis=1)))
            sampled_xy[inside] = np.einsum(
                'ni,nij->nj', weights, src[triangulation.simplices[ids]])
        # Output pixel samples the position in the input plane: (x,y)+D(x,y).
        values = ndimage.map_coordinates(
            source_mask.astype(np.uint8),
            [sampled_xy[:, 1] - 1, sampled_xy[:, 0] - 1],
            order=0, mode="constant", cval=0, prefilter=False)
        warped[top:bottom, x_start:x_end] = values.reshape(bottom - top, x_end - x_start) > 0
    return warped


def add_outline_to_thumbnail(thumb: Path, mask: np.ndarray, destination: Path) -> None:
    """Preserve original histology; draw ONLY the Su LC perimeter."""
    with Image.open(thumb) as opened:
        rgb = np.array(opened.convert("RGB"))
    if rgb.shape[:2] != mask.shape:
        raise ValueError(f"Thumbnail {rgb.shape[:2]} and LC {mask.shape} size mismatch")
    perimeter = mask & ~ndimage.binary_erosion(mask, structure=np.ones((3, 3), bool))
    if OUTLINE_THICKNESS > 1:
        perimeter = ndimage.binary_dilation(perimeter, iterations=OUTLINE_THICKNESS - 1)
    rgb[perimeter] = OUTLINE_COLOR
    Image.fromarray(rgb, "RGB").save(destination)


def find_thumbnail(mouse_folder: Path, name: str) -> Path:
    folder = mouse_folder / "thumbnails"
    preferred = folder / (name + "-thumb.png")
    if preferred.exists():
        return preferred
    matches = sorted(folder.glob(name + "*.png"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one thumbnail for {name}; found {len(matches)} in {folder}")
    return matches[0]


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    mouse_folder = PROJECT_FOLDER / "rawData" / MOUSE_ID
    info = mouse_folder / f"{MOUSE_ID}-info.xml"
    original_json = mouse_folder / f"{MOUSE_ID}-visualign.json"
    mesh_file = PROJECT_FOLDER / MESH_RELATIVE_PATH
    output = PROJECT_FOLDER / "LC_analysis" / MOUSE_ID / "SuLC_only"
    for path, desc in ((info, "Mouse info XML"),
                       (original_json, "Original VisuAlign JSON"),
                       (mesh_file, "Su LC core mesh")):
        _require(path, desc)

    sections = load_sections(info)
    records = load_alignment_records(original_json)
    q, faces = load_bilateral_quicknii_mesh(mesh_file)
    output.mkdir(parents=True, exist_ok=True)
    print(f"Mouse: {MOUSE_ID} | sections: {len(sections)}")
    print(f"Su mesh: {len(q)} bilateral vertices, {len(faces)} triangles")
    print(f"Output: {output}")

    fields = ["MouseID", "SliceName", "SliceNumber", "Valid", "HasSuLC",
              "SuPixelsInPlane", "SuPixelsAfterOriginalWarp", "MeshSegments",
              "ImageFile", "Status"]
    rows = []
    for idx, sec in enumerate(sections, start=1):
        name = sec["name"]
        row = dict(MouseID=MOUSE_ID, SliceName=name, SliceNumber=sec["number"],
                   Valid=sec["valid"], HasSuLC=False, SuPixelsInPlane=0,
                   SuPixelsAfterOriginalWarp=0, MeshSegments=0,
                   ImageFile="", Status="")
        print(f"[{idx:3d}/{len(sections)}] {name}: ", end="", flush=True)
        if SKIP_INVALID_SLICES and not sec["valid"]:
            row["Status"] = "Skipped invalid section"
        else:
            try:
                thumb = find_thumbnail(mouse_folder, name)
                with Image.open(thumb) as img:
                    w, h = img.size
                if name not in records:
                    raise ValueError("No exact match in original VisuAlign JSON")
                rec = records[name]
                anchoring = np.asarray(rec.get("anchoring", []), dtype=float).ravel()
                if anchoring.shape != (9,) or not np.isfinite(anchoring).all():
                    raise ValueError("Original QuickNII anchoring must have 9 finite values")

                segments = intersect_mesh_with_plane(q, faces, anchoring, w, h)
                row["MeshSegments"] = len(segments)
                plane_mask = rasterize_section(segments, h, w)
                n_plane = int(np.count_nonzero(plane_mask))
                row["SuPixelsInPlane"] = n_plane
                if n_plane < MIN_SU_PIXELS:
                    row["Status"] = ("Surface intersects but no rasterized pixels" if len(segments)
                                     else "No Su LC intersection")
                else:
                    # Detect 3-D intersection FIRST; VisuAlign is an in-plane warp.
                    # Same selection principle as find_Su_LC_slices.m.
                    row["HasSuLC"] = True
                    markers = np.asarray(rec.get("markers", []), dtype=float)
                    aligned_mask = warp_mask_with_original_markers(plane_mask, markers)
                    n_aligned = int(np.count_nonzero(aligned_mask))
                    row["SuPixelsAfterOriginalWarp"] = n_aligned
                    if not n_aligned:
                        row["Status"] = "ERROR: Su LC vanished after original warp"
                    else:
                        image_name = f"{name}_SuLCOverlay.png"
                        add_outline_to_thumbnail(thumb, aligned_mask, output / image_name)
                        row["ImageFile"] = image_name
                        row["Status"] = "Saved"
            except Exception as exc:
                row["Status"] = f"ERROR: {exc}"
        rows.append(row)
        print(f"{row['Status']} | plane={row['SuPixelsInPlane']} "
              f"aligned={row['SuPixelsAfterOriginalWarp']}")

    all_path = output / f"{MOUSE_ID}_SuLC_all_slices.csv"
    positive_path = output / f"{MOUSE_ID}_SuLC_positive_slices.csv"
    write_csv(all_path, rows, fields)
    write_csv(positive_path, [r for r in rows if r["HasSuLC"]], fields)
    saved = sum(r["Status"] == "Saved" for r in rows)
    detected = sum(r["HasSuLC"] for r in rows)
    errors = sum(r["Status"].startswith("ERROR:") for r in rows)
    print(f"\nComplete: detected {detected}, saved {saved} PNGs, errors {errors}.")
    print(f"All-slice CSV: {all_path}")
    print(f"Positive CSV:  {positive_path}")
    if errors:
        print("Review error rows in the all-slices CSV (no input files were modified).")


if __name__ == "__main__":
    main()

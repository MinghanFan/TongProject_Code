%% Inspect Su et al. LC mesh in relation to Allen CCFv3 LC
% Diagnostic only: checks spatial orientation; does not edit atlas or JSON.
clear; clc; close all;

projectFolder = '/Users/minghanfan/Desktop/TongProject/TongProject_Code';
suFolder = fullfile(projectFolder, 'atlas', 'Su_LC');
meshFile = fullfile(suFolder, 'core_mesh.obj');
csvFile = fullfile(suFolder, 'lc_points_with_features.csv');
atlasFile = fullfile(projectFolder, 'atlas', 'annotationVolume.mat');

assert(isfile(meshFile), 'Cannot find mesh: %s', meshFile);
assert(isfile(csvFile), 'Cannot find CSV: %s', csvFile);
assert(isfile(atlasFile), 'Cannot find atlas: %s', atlasFile);

%% CSV: inspect core membership, not just all cropped points
T = readtable(csvFile);
corePoints = T(T.in_core_mesh ~= 0, :);
fprintf('Cropped point calls: %d\n', height(T));
fprintf('Points inside core mesh: %d\n', height(corePoints));
if ~isempty(corePoints)
    fprintf('Core-points range (um), x/y/z:\n');
    disp([min(corePoints{:, {'x','y','z'}}, [], 1); ...
          max(corePoints{:, {'x','y','z'}}, [], 1)]);
end

%% Su core mesh (one midline-folded hemisphere)
% readSurfaceMesh supports Wavefront OBJ since MATLAB R2024a, when available.
if exist('readSurfaceMesh', 'file')
    mesh = readSurfaceMesh(meshFile);
    V = double(mesh.Vertices);
    F = double(mesh.Faces);
else
    % Fallback OBJ parser if readSurfaceMesh is not installed.
    [V, F] = readSimpleOBJ(meshFile);
end
fprintf('Core mesh: %d vertices, %d triangular faces\n', size(V,1), size(F,1));
fprintf('Mesh vertex range (um), x/y/z:\n');
disp([min(V, [], 1); max(V, [], 1)]);

%% Reconstruct the contralateral copy across the authors' x=5700 um midline
Vother = V;
Vother(:,1) = 11400 - Vother(:,1);
Vboth = [V; Vother];
Fboth = [F; F + size(V,1)];

%% Existing CCFv3 LC from your project's atlas preprocessing
S = load(atlasFile);
assert(isfield(S, 'annotationVolume'), 'MAT file lacks annotationVolume.');
A = permute(S.annotationVolume, [3 1 2]);
A = flip(A, 2);
A = flip(A, 3);

[ii,jj,kk] = ind2sub(size(A), find(A == 147));
assert(~isempty(ii), 'No old CCFv3 LC voxels (ID 147) found.');

% Candidate physical-coordinate transform corresponding to the atlas flips.
% Verify visually against Su's bilateral mesh; do not treat as final until checked.
oldXYZ_um = 25 * [double(ii)-1, ...
                  double(size(A,2))-double(jj), ...
                  double(size(A,3))-double(kk)];

fprintf('Original CCFv3 LC: %d voxels\n', size(oldXYZ_um,1));
fprintf('Old LC candidate range (um), x/y/z:\n');
disp([min(oldXYZ_um, [], 1); max(oldXYZ_um, [], 1)]);

%% View alignment in 3D
fig = figure('Color','w', 'Position', [100 100 1100 760]);
ax = axes(fig);
patch(ax, 'Faces', Fboth, 'Vertices', Vboth, ...
    'FaceColor', [0.0 0.75 0.82], 'FaceAlpha', 0.22, ...
    'EdgeColor', 'none');
hold(ax, 'on');
scatter3(ax, oldXYZ_um(:,1), oldXYZ_um(:,2), oldXYZ_um(:,3), ...
    9, [0.9 0.1 0.1], 'filled');
axis(ax, 'equal'); grid(ax, 'on'); view(ax, 3);
xlabel(ax, 'Su x / ML (um)');
ylabel(ax, 'Su y / AP (um)');
zlabel(ax, 'Su z / DV (um)');
title(ax, 'Spatial check: Su bilateral core (cyan) vs CCFv3 LC 147 (red)');
legend(ax, {'Su core meshes', 'CCFv3 LC voxels'}, 'Location', 'best');

outFile = fullfile(suFolder, 'su_vs_ccfv3_LC_spatial_check.png');
exportgraphics(fig, outFile, 'Resolution', 180);
fprintf('\nSaved diagnostic figure to:\n%s\n', outFile);
fprintf('CHECK: Do the red LC voxels appear near or inside the cyan meshes?\n');

%% Local helper: basic OBJ reader (supports triangular faces)
function [V,F] = readSimpleOBJ(path)
    fid = fopen(path, 'r');
    assert(fid ~= -1, 'Could not open %s', path);
    lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = lines{1};
    vl = lines(startsWith(lines, 'v '));
    fl = lines(startsWith(lines, 'f '));
    V = zeros(numel(vl), 3);
    F = zeros(numel(fl), 3);
    for q = 1:numel(vl)
        nums = sscanf(vl{q}(3:end), '%f');
        V(q,:) = nums(1:3)';
    end
    for q = 1:numel(fl)
        parts = regexp(strtrim(fl{q}(3:end)), '\s+', 'split');
        assert(numel(parts)==3, 'This OBJ contains non-triangular faces.');
        for j = 1:3
            first = regexp(parts{j}, '^-?\d+', 'match', 'once');
            F(q,j) = str2double(first);
        end
    end
    assert(all(F(:)>=1), 'OBJ uses unsupported negative/invalid face indices.');
end

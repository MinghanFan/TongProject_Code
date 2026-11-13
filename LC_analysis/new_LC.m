%% Set 3: Manual VisuAlign CCFv3 atlas + Su et al. LC core mesh
% Generates: histology.png, atlas.png, overlay.png.
% Does NOT modify source JSON, XML, atlas, or displacement CSV files.
% Uses an exact plane/surface intersection of the published LC core mesh,
% NOT a dilated LC point cloud or an alpha-shape approximation.
clear; clc; close all;

%% USER SETTINGS
projectFolder = '/Users/minghanfan/Desktop/TongProject/TongProject_Code';
mouseID = 'CC1A';
sliceName = 'CC1A_172_B2_1';

addpath(genpath(fullfile(projectFolder, 'brainAlignment')));
mouseFolder = fullfile(projectFolder, 'rawData', mouseID);
infoFile = fullfile(mouseFolder, [mouseID '-info.xml']);
manualJson = fullfile(mouseFolder, [mouseID '-visualign-manual.json']);
atlasFile = fullfile(projectFolder, 'atlas', 'annotationVolume.mat');
meshFile = fullfile(projectFolder, 'atlas', 'Su_LC', 'core_mesh.obj');

outDir = fullfile(projectFolder, 'LC_analysis', mouseID, ...
    'comparison', sliceName, '03_Manual_VisuAlign_SuLC');
if ~isfolder(outDir); mkdir(outDir); end

lcID = 147;            % Original CCFv3 LC label
voxelSizeUm = 25;      % Your CCFv3 annotation volume resolution
midlineUm = 5700;     % Su mesh was constructed with both sides folded at x=5700 um
resizeFactor = 0.2;   % Reuse the value tested for the manual alignment

%% 1. LOAD THE SAME MANUALLY REFINED REGISTRATION AS SET 2
assert(isfile(infoFile), 'Cannot find %s', infoFile);
assert(isfile(manualJson), 'Cannot find %s', manualJson);
assert(isfile(atlasFile), 'Cannot find %s', atlasFile);
assert(isfile(meshFile), 'Cannot find %s', meshFile);

sManual = Slice(sliceName, infoFile);
assert(~isempty(sManual.thumbnail), 'Thumbnail is empty.');
record = findSliceRecord(manualJson, sliceName);
a = double(record.anchoring(:)');
assert(numel(a)==9, 'Expected 9 QuickNII anchoring parameters.');
sManual.o = a(1:3);
sManual.u = a(4:6);
sManual.v = a(7:9);
sManual.markers = double(record.markers);
assert(size(sManual.markers,2)==4 && ~isempty(sManual.markers), ...
    'Expected nonempty Nx4 markers in the manual JSON.');
D = sManual.getDispField(resizeFactor);
sManual.dispField = imresize(D, size(sManual.thumbnail,[1 2]), 'bilinear');

%% 2. LOAD CCFv3 WITH THE SAME ORIENTATION AS YOUR PREVIOUS MATLAB SCRIPTS
S = load(atlasFile, 'annotationVolume');
A = permute(S.annotationVolume, [3 1 2]);
A = flip(A, 2);
A = flip(A, 3);

% Original atlas labels, using exactly the same QuickNII and VisuAlign transforms.
manualPlane = sManual.getObliqueSlice(A);
manualAtlas = imwarp(manualPlane, sManual.dispField, 'nearest');
[H,W,~] = size(sManual.thumbnail);
assert(isequal(size(manualAtlas),[H W]), 'Atlas size differs from thumbnail.');
oldLC = (manualAtlas == lcID);

%% 3. READ THE AUTHOR'S PUBLISHED LC CORE MESH (WAVEFRONT OBJ)
if exist('readSurfaceMesh','file')
    mesh = readSurfaceMesh(meshFile);
    V = double(mesh.Vertices);        % Su coordinates, in micrometers
    F = double(mesh.Faces);           % Triangular faces
else
    [V,F] = readSimpleOBJ(meshFile);
end
assert(size(V,2)==3 && size(F,2)==3, 'Expected triangle mesh.');
assert(all(F(:)>=1) && all(F(:)<=size(V,1)), 'Invalid face indices.');

% The paper's reconstruction folds both sides onto x <= 5700 um.
% Reconstruct a separate opposite-hemisphere copy in the same CCF space.
Vright = V;
Vright(:,1) = 2*midlineUm - Vright(:,1);
Vbilateral = [V; Vright];
Fbilateral = [F; F+size(V,1)];

%% 4. CONVERT SU COORDINATES TO QUICKNII ATLAS COORDINATES
% The mapping below is derived from the 3D check performed previously:
%  Su physical XYZ (um) = 25*[i-1, size(A,2)-j, size(A,3)-k]
% for MATLAB atlas voxel indices (i,j,k).
%
% IMPORTANT: Slice.getObliqueSlice indexes atlas samples as
% b = round(o + a*u + b*v) + [0 1 1].
% Thus its continuous QuickNII coordinate values for the mesh are:
Q = [1 + Vbilateral(:,1)/voxelSizeUm, ...
     size(A,2)-1 - Vbilateral(:,2)/voxelSizeUm, ...
     size(A,3)-1 - Vbilateral(:,3)/voxelSizeUm];

%% 5. INTERSECT THE 3D TRIANGLE MESH WITH THE ACTUAL QUICKNII SECTION PLANE
% This gives the Su LC cross-section in thumbnail space, before VisuAlign.
segments = intersectMeshWithQuickNiiPlane(Q,Fbilateral, ...
    sManual.o,sManual.u,sManual.v,W,H);
fprintf('Mesh/QuickNII plane intersection: %d segments\n',size(segments,1));

suPlaneMask = rasterizeClosedMeshSections(segments,H,W);
% Warp using the EXACT SAME manual VisuAlign displacement field as the atlas.
suAligned = imwarp(suPlaneMask,sManual.dispField,'nearest');
assert(isequal(size(suAligned),[H W]), 'Su mask size differs from thumbnail.');

fprintf('Original CCFv3 LC pixels (manual alignment): %d\n',nnz(oldLC));
fprintf('Su LC core pixels (manual alignment):       %d\n',nnz(suAligned));
if any(oldLC(:))
    fprintf('Original LC pixels also in Su LC core: %.1f%%\n', ...
        100*nnz(oldLC & suAligned)/nnz(oldLC));
end
if ~any(suAligned(:))
    warning(['Su LC has no intersection with this QuickNII section. ' ...
             'Confirm the plane and Su coordinate mapping before interpreting the output.']);
end

%% 6. SAVE THE THREE REQUESTED IMAGES
rgb = sManual.thumbnail;
if size(rgb,3)==1; rgb = repmat(rgb,[1 1 3]); end
if ~isa(rgb,'uint8'); rgb = im2uint8(rgb); end
imwrite(rgb,fullfile(outDir,'histology.png'));

borders = Slice.annotationBorders(Slice.rescaleAnnotation(manualAtlas));
oldEdge = bwperim(oldLC,8);
suEdge = bwperim(suAligned,8);

% Atlas only: white CCFv3 boundaries, red original LC, cyan Su core.
atlasView = uint8(25*ones(H,W,3));
atlasView = blendPixels(atlasView,suAligned,[0 205 220],0.27);
atlasView = blendPixels(atlasView,borders,[220 220 220],1.00);
atlasView = blendPixels(atlasView,oldEdge,[255 55 55],1.00);
atlasView = blendPixels(atlasView,suEdge,[0 255 255],1.00);
imwrite(atlasView,fullfile(outDir,'atlas.png'));

% Histology overlay: identical manual atlas, with Su mesh shown explicitly.
overlay = blendPixels(rgb,suAligned,[0 205 220],0.25);
overlay = blendPixels(overlay,borders,[255 255 255],0.50);
overlay = blendPixels(overlay,oldEdge,[255 45 45],1.00);
overlay = blendPixels(overlay,suEdge,[0 255 255],1.00);
imwrite(overlay,fullfile(outDir,'overlay.png'));

fprintf('\nSaved Set 3 to:\n%s\n',outDir);
fprintf('Colors: red = CCFv3 LC 147; cyan = Su LC core; white = CCFv3 borders.\n');
fprintf('No original registration/data files were overwritten.\n');

%% LOCAL HELPER FUNCTIONS
function record = findSliceRecord(jsonFile, sliceName)
    J = jsondecode(strrep(fileread(jsonFile),'\','/'));
    assert(isfield(J,'slices'), 'No slices in JSON.');
    files = string({J.slices.filename});
    hit = find(contains(files,sliceName));
    assert(numel(hit)==1, 'Expected one JSON record for %s; found %d.', ...
        sliceName,numel(hit));
    record = J.slices(hit);
end

function seg = intersectMeshWithQuickNiiPlane(Q,F,o,u,v,W,H)
    % Slice the triangular surface against p = o + alpha*u + beta*v.
    % Return 2D line segments [x1,y1,x2,y2] in thumbnail coordinates.
    n = cross(u,v);
    assert(norm(n)>0, 'QuickNII plane u,v vectors are degenerate.');
    n = n/norm(n);
    d = (Q - double(o(:)'))*n(:);
    triD = reshape(d(F(:)),size(F));
    chosen = (min(triD,[],2)<=0) & (max(triD,[],2)>0);
    Fc = F(chosen,:);
    if isempty(Fc); seg=zeros(0,4); return; end
    entries = zeros(0,3);
    basis = [u(:),v(:)];
    for e=1:3
        e1=e; e2=mod(e,3)+1;
        ia=Fc(:,e1); ib=Fc(:,e2);
        da=d(ia); db=d(ib);
        crossEdge = (da<=0 & db>0) | (db<=0 & da>0);
        if ~any(crossEdge); continue; end
        ta=find(crossEdge);
        aa=ia(crossEdge); bb=ib(crossEdge);
        f=da(crossEdge)./(da(crossEdge)-db(crossEdge));
        pts=Q(aa,:)+f.*(Q(bb,:)-Q(aa,:));
        coeff=basis\(pts-double(o(:)')).';
        xy=[coeff(1,:)'*W+1, coeff(2,:)'*H+1];
        entries=[entries; [ta,xy]]; %#ok<AGROW>
    end
    if isempty(entries); seg=zeros(0,4); return; end
    entries=sortrows(entries,1);
    counts=accumarray(entries(:,1),1,[size(Fc,1),1]);
    keep=(counts(entries(:,1))==2);
    if any(~keep)
        warning('Skipped %d ambiguous surface/plane intersection entries.',nnz(~keep));
    end
    entries=entries(keep,:);
    if isempty(entries); seg=zeros(0,4); return; end
    assert(mod(size(entries,1),2)==0,'Uneven triangle/plane crossings.');
    seg=reshape(entries(:,2:3).',4,[]).';
    nonzero=hypot(seg(:,1)-seg(:,3),seg(:,2)-seg(:,4))>1e-8;
    seg=seg(nonzero,:);
end

function mask = rasterizeClosedMeshSections(seg,H,W)
    % Even/odd rule for all closed loops in the planar surface intersection.
    % Handles multiple disconnected components without reconstructing loops.
    mask=false(H,W);
    if isempty(seg); return; end
    lowRow=max(1,ceil(min(min(seg(:,[2 4])))));
    highRow=min(H,floor(max(max(seg(:,[2 4])))));
    for r=lowRow:highRow
        y0=seg(:,2); y1=seg(:,4);
        crossRow=(y0<=r & y1>r) | (y1<=r & y0>r);
        if ~any(crossRow); continue; end
        s=seg(crossRow,:);
        xc=s(:,1)+(r-s(:,2)).*(s(:,3)-s(:,1))./(s(:,4)-s(:,2));
        xc=sort(xc);
        if mod(numel(xc),2)~=0
            warning('Odd surface crossings in row %d; row omitted.',r);
            continue;
        end
        for k=1:2:numel(xc)
            xleft=max(1,ceil(xc(k)));
            xright=min(W,floor(xc(k+1)));
            if xright>=xleft
                mask(r,xleft:xright)=true;
            end
        end
    end
end

function img = blendPixels(img,mask,color,alpha)
    img=double(img);
    for c=1:3
        ch=img(:,:,c);
        ch(mask)=(1-alpha)*ch(mask)+alpha*color(c);
        img(:,:,c)=ch;
    end
    img=uint8(round(img));
end

function [V,F] = readSimpleOBJ(path)
    % Simple OBJ reader, for triangle-mesh files if readSurfaceMesh unavailable.
    fid=fopen(path,'r');
    assert(fid~=-1, 'Unable to open OBJ file.');
    raw=textscan(fid,'%s','Delimiter','\n','Whitespace','');
    fclose(fid);
    lines=raw{1};
    vlines=lines(startsWith(lines,'v '));
    flines=lines(startsWith(lines,'f '));
    V=zeros(numel(vlines),3);
    F=zeros(numel(flines),3);
    for k=1:numel(vlines)
        n=sscanf(vlines{k}(3:end),'%f');
        V(k,:)=n(1:3)';
    end
    for k=1:numel(flines)
        parts=regexp(strtrim(flines{k}(3:end)),'\s+','split');
        assert(numel(parts)==3,'OBJ has a nontriangular face.');
        for j=1:3
            token=regexp(parts{j},'^\d+','match','once');
            assert(~isempty(token),'OBJ face contains a negative/invalid index.');
            F(k,j)=str2double(token);
        end
    end
end

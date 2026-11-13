%% FIND Su et al. LC CORE IN ALL SECTIONS OF ONE MOUSE
% Run this script FIRST. Change only mouseID to scan a different mouse.
%
% Input:  CC1A-info.xml, CC1A-visualign-manual.json, Su core_mesh.obj
% Output: CC1A_SuLC_all_slices.csv, CC1A_SuLC_positive_slices.csv
%
% IMPORTANT: Su LC is detected in the MANUAL QuickNII 3D section plane,
% BEFORE the 2D VisuAlign displacement. In-plane warping cannot change
% which 3D anatomical plane is sampled. The export script applies the
% manual displacement to the actual overlays.
clear; clc;

%% USER SETTINGS
mouseID = 'CC1A';                 % Change this for another mouse
projectFolder = '/Users/minghanfan/Desktop/TongProject/TongProject_Code';
scanInvalidSlices = false;       % True = include invalid histology sections
minSuPixels = 1;                 % 1 = count any rasterized LC core pixel

addpath(genpath(fullfile(projectFolder,'brainAlignment')));
addpath(fileparts(mfilename('fullpath'))); % Finds su_lc_tools.m
mouseFolder = fullfile(projectFolder,'rawData',mouseID);
infoFile = fullfile(mouseFolder,[mouseID '-info.xml']);
manualJson = fullfile(mouseFolder,[mouseID '-visualign-manual.json']);
meshFile = fullfile(projectFolder,'atlas','Su_LC','core_mesh.obj');
atlasFile = fullfile(projectFolder,'atlas','annotationVolume.mat');
outFolder = fullfile(projectFolder,'LC_analysis',mouseID,'SuLC_batch');

voxelSizeUm = 25;
midlineUm = 5700;

assert(isfile(infoFile),'Info XML not found: %s',infoFile);
assert(isfile(manualJson),'Manual JSON not found: %s',manualJson);
assert(isfile(meshFile),'Su mesh not found: %s',meshFile);
assert(isfile(atlasFile),'Atlas not found: %s',atlasFile);
if ~isfolder(outFolder); mkdir(outFolder); end

%% LOAD METADATA AND SU MESH JUST ONCE
X = readstruct(infoFile);
sections = X.slices;
n = numel(sections);
names = string({sections.name})';
numbers = double([sections.number])';
valid = logical([sections.valid])';

records = su_lc_tools('loadRecords',manualJson);
atlasInfo = whos('-file',atlasFile,'annotationVolume');
assert(~isempty(atlasInfo),'annotationVolume not found in MAT file.');
atlasSize = atlasInfo.size([3 1 2]);   % permute([3 1 2]); flips preserve size
[Q,F] = su_lc_tools('loadMesh',meshFile,atlasSize,voxelSizeUm,midlineUm);
fprintf('Loaded Su bilateral LC mesh: %d vertices, %d triangles.\n', ...
    size(Q,1),size(F,1));

%% SCAN EACH SECTION WITHOUT READING FULL ATLAS OR DISPLACEMENT CSVS
hasSuLC = false(n,1);
suPixels = zeros(n,1);
meshSegments = zeros(n,1);
thumbWidth = zeros(n,1);
thumbHeight = zeros(n,1);
status = strings(n,1);

for i = 1:n
    name = char(names(i));
    fprintf('[%d/%d] %-24s ',i,n,name);
    if ~valid(i) && ~scanInvalidSlices
        status(i) = "Skipped: invalid slice";
        fprintf('%s\n',status(i));
        continue;
    end
    try
        % Read only the thumbnail dimensions (not the entire image).
        thumbnails = dir(fullfile(mouseFolder,'thumbnails',[name '*.png']));
        assert(~isempty(thumbnails),'Missing thumbnail PNG.');
        thumbFile = fullfile(thumbnails(1).folder,thumbnails(1).name);
        im = imfinfo(thumbFile);
        W = double(im.Width);
        H = double(im.Height);
        thumbWidth(i) = W;
        thumbHeight(i) = H;

        rec = su_lc_tools('findRecord',records,name);
        a = double(rec.anchoring(:)');
        assert(numel(a)==9,'Expected nine QuickNII anchoring values.');
        [M,nseg] = su_lc_tools('sliceMask',Q,F,...
            a(1:3),a(4:6),a(7:9),W,H);
        suPixels(i) = nnz(M);
        meshSegments(i) = nseg;
        hasSuLC(i) = suPixels(i)>=minSuPixels;
        if hasSuLC(i)
            status(i) = "Su LC found";
        elseif nseg>0
            status(i) = "Surface intersects; no filled pixels";
        else
            status(i) = "No Su LC";
        end
        fprintf('%s | pixels=%d | segments=%d\n', ...
            status(i),suPixels(i),nseg);
    catch ME
        status(i) = "ERROR: " + string(ME.message);
        warning('Scan failed for %s: %s',name,ME.message);
    end
end

%% WRITE FULL REFERENCE AND POSITIVE-SLICES CSV
R = table(repmat(string(mouseID),n,1),names,numbers,valid,hasSuLC,...
    suPixels,meshSegments,thumbWidth,thumbHeight,status,...
    'VariableNames',{'MouseID','SliceName','SliceNumber','Valid', ...
    'HasSuLC','SuPixelsInPlane','MeshSegments','Width','Height','Status'});
allCsv = fullfile(outFolder,[mouseID '_SuLC_all_slices.csv']);
hitsCsv = fullfile(outFolder,[mouseID '_SuLC_positive_slices.csv']);
writetable(R,allCsv);
writetable(R(R.HasSuLC,:),hitsCsv);

fprintf('\nDONE: %d / %d sections contain Su LC (>= %d pixels).\n', ...
    nnz(R.HasSuLC),n,minSuPixels);
fprintf('All sections:      %s\n',allCsv);
fprintf('Positive sections: %s\n',hitsCsv);
fprintf('Errors: %d; inspect Status column in all-slices CSV.\n', ...
    nnz(startsWith(status,"ERROR:")));
fprintf('No input files or existing displacement CSVs were changed.\n');

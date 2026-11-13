%% EXPORT ALL THREE ALIGNMENT SETS FOR Su-LC-POSITIVE SECTIONS
% Run find_Su_LC_slices.m first for the same mouseID.
% This script reads the positive-slices CSV and creates NINE PNGs per slice:
%   01_Original/             histology.png, atlas.png, overlay.png
%   02_Manual_VisuAlign/     histology.png, atlas.png, overlay.png
%   03_Manual_VisuAlign_SuLC/ histology.png, atlas.png, overlay.png
%
% Your original JSON, manual JSON, XML and displacement CSVs are preserved.
clear; clc; close all;

%% USER SETTINGS
mouseID = 'CC1A';                 % Change this for another mouse
projectFolder = '/Users/minghanfan/Desktop/TongProject/TongProject_Code';
resizeFactor = 0.2;              % Same tested VisuAlign marker sampling
lcID = 147;                      % Original Allen CCFv3 LC region ID

% Normally export EVERY Su-positive section in the reference CSV.
% Optional: restrict to selected names, for example:
% onlyTheseSlices = ["CC1A_160_A4_1", "CC1A_163_A5_1"];
onlyTheseSlices = strings(0,1);
maxSlices = inf;                 % Set to 1 for a quick trial, inf for all

voxelSizeUm = 25;
midlineUm = 5700;

addpath(genpath(fullfile(projectFolder,'brainAlignment')));
addpath(fileparts(mfilename('fullpath'))); % Finds su_lc_tools.m
mouseFolder = fullfile(projectFolder,'rawData',mouseID);
infoFile = fullfile(mouseFolder,[mouseID '-info.xml']);
originalJson = fullfile(mouseFolder,[mouseID '-visualign.json']);
manualJson = fullfile(mouseFolder,[mouseID '-visualign-manual.json']);
atlasFile = fullfile(projectFolder,'atlas','annotationVolume.mat');
meshFile = fullfile(projectFolder,'atlas','Su_LC','core_mesh.obj');
rootOut = fullfile(projectFolder,'LC_analysis',mouseID,'SuLC_batch');
refCsv = fullfile(rootOut,[mouseID '_SuLC_positive_slices.csv']);

assert(isfile(refCsv),['Reference CSV not found: %s\n' ...
    'Run find_Su_LC_slices.m for this mouse first.'],refCsv);
assert(isfile(infoFile),'Info XML not found: %s',infoFile);
assert(isfile(originalJson),'Original JSON not found: %s',originalJson);
assert(isfile(manualJson),'Manual JSON not found: %s',manualJson);
assert(isfile(atlasFile),'Atlas not found: %s',atlasFile);
assert(isfile(meshFile),'Su core mesh not found: %s',meshFile);

%% LOAD CSV AND SELECT SECTIONS
R = readtable(refCsv,'TextType','string');
required = {'SliceName','HasSuLC','SuPixelsInPlane'};
assert(all(ismember(required,R.Properties.VariableNames)), ...
    'Reference CSV is missing required columns.');
if islogical(R.HasSuLC)
    hit = R.HasSuLC;
elseif isnumeric(R.HasSuLC)
    hit = R.HasSuLC~=0;
else
    hit = strcmpi(string(R.HasSuLC),"true") | string(R.HasSuLC)=="1";
end
R = R(hit,:);
if ~isempty(onlyTheseSlices)
    R = R(ismember(R.SliceName,onlyTheseSlices),:);
end
R = R(1:min(height(R),maxSlices),:);
assert(~isempty(R),'No Su-LC-positive slices selected from %s.',refCsv);

%% LOAD EACH RESOURCE ONCE
recordsOriginal = su_lc_tools('loadRecords',originalJson);
recordsManual = su_lc_tools('loadRecords',manualJson);
S = load(atlasFile,'annotationVolume');
A = permute(S.annotationVolume,[3 1 2]);
A = flip(A,2);
A = flip(A,3);
[Q,F] = su_lc_tools('loadMesh',meshFile,size(A),voxelSizeUm,midlineUm);

%% EXPORT SELECTED SECTIONS
n = height(R);
originalLCPixels = zeros(n,1);
manualLCPixels = zeros(n,1);
suPlanePixels = zeros(n,1);
suAlignedPixels = zeros(n,1);
exportStatus = strings(n,1);
outputDir = strings(n,1);

for i = 1:n
    sliceName = char(R.SliceName(i));
    fprintf('\n[%d/%d] EXPORT %s\n',i,n,sliceName);
    try
        % Slice supplies histology + dataset's ORIGINAL displacement CSV.
        sBase = Slice(sliceName,infoFile);
        assert(~isempty(sBase.thumbnail),'Thumbnail is missing.');
        [H,W,~] = size(sBase.thumbnail);

        % Explicitly select original and manual JSONs (Slice.m only
        % auto-detects files ending in '-visualign.json').
        recOrig = su_lc_tools('findRecord',recordsOriginal,sliceName);
        recManual = su_lc_tools('findRecord',recordsManual,sliceName);
        sOriginal = su_lc_tools('applyRecord',sBase,recOrig);
        sManual = su_lc_tools('applyRecord',sBase,recManual);

        if isempty(sBase.dispField)
            assert(~isempty(sOriginal.markers), ...
                'No original CSV displacement or original markers.');
            d0 = sOriginal.getDispField(resizeFactor);
            sOriginal.dispField = imresize(d0,[H W],'bilinear');
        else
            % Retain the dataset's original displacement field for Set 1.
            sOriginal.dispField = sBase.dispField;
        end
        assert(~isempty(sManual.markers), ...
            'Manual JSON has no VisuAlign markers for this slice.');
        dm = sManual.getDispField(resizeFactor);
        sManual.dispField = imresize(dm,[H W],'bilinear');

        % Same CCF volume, sampled with each version's QuickNII vectors.
        originalPlane = sOriginal.getObliqueSlice(A);
        manualPlane = sManual.getObliqueSlice(A);
        originalAtlas = imwarp(originalPlane,sOriginal.dispField,'nearest');
        manualAtlas = imwarp(manualPlane,sManual.dispField,'nearest');
        assert(isequal(size(originalAtlas),[H W]) && ...
               isequal(size(manualAtlas),[H W]), ...
               'Unexpected output size after VisuAlign warp.');

        % Su LC: intersect the real core mesh with the MANUAL plane,
        % then warp the mask with the SAME manual displacement as Set 2.
        [suPlane,~] = su_lc_tools('sliceMask',Q,F, ...
            sManual.o,sManual.u,sManual.v,W,H);
        suAligned = imwarp(suPlane,sManual.dispField,'nearest');
        assert(isequal(size(suAligned),[H W]),'Su mask size differs.');
        if nnz(suPlane)~=R.SuPixelsInPlane(i)
            warning(['Su pixels differ from scan CSV for %s (%d now vs %d before). ' ...
                'If you edited the manual JSON, re-run find_Su_LC_slices.m.'], ...
                sliceName,nnz(suPlane),R.SuPixelsInPlane(i));
        end

        folder = fullfile(rootOut,'final_outputs',sliceName);
        d1 = fullfile(folder,'01_Original');
        d2 = fullfile(folder,'02_Manual_VisuAlign');
        d3 = fullfile(folder,'03_Manual_VisuAlign_SuLC');
        su_lc_tools('saveImages',sBase.thumbnail,originalAtlas,lcID,d1);
        su_lc_tools('saveImages',sBase.thumbnail,manualAtlas,lcID,d2);
        su_lc_tools('saveImages',sBase.thumbnail,manualAtlas,lcID,d3,suAligned);

        originalLCPixels(i) = nnz(originalAtlas==lcID);
        manualLCPixels(i) = nnz(manualAtlas==lcID);
        suPlanePixels(i) = nnz(suPlane);
        suAlignedPixels(i) = nnz(suAligned);
        outputDir(i) = string(folder);
        exportStatus(i) = "Saved";
        if ~any(suAligned(:))
            exportStatus(i) = "Saved, but Su LC became empty after warp";
            warning('Su LC empty after manual warp: %s',sliceName);
        end
        fprintf('Original LC=%d, Manual LC=%d, Su plane=%d, Su aligned=%d\n', ...
            originalLCPixels(i),manualLCPixels(i), ...
            suPlanePixels(i),suAlignedPixels(i));
    catch ME
        exportStatus(i) = "ERROR: " + string(ME.message);
        warning('Export failed for %s: %s',sliceName,ME.message);
    end
end

%% SAVE AN AUDIT LOG OF WHAT WAS EXPORTED
L = table(R.SliceName,originalLCPixels,manualLCPixels,...
    suPlanePixels,suAlignedPixels,exportStatus,outputDir,...
    'VariableNames',{'SliceName','OriginalCCFLCPixels','ManualCCFLCPixels', ...
    'SuPixelsInPlane','SuPixelsAfterManualWarp','ExportStatus','OutputFolder'});
logCsv = fullfile(rootOut,[mouseID '_SuLC_export_summary.csv']);
writetable(L,logCsv);
fprintf('\nDONE. %d / %d selected sections saved.\n', ...
    nnz(startsWith(exportStatus,"Saved")),n);
fprintf('Images:  %s\n',fullfile(rootOut,'final_outputs'));
fprintf('Summary: %s\n',logCsv);
fprintf('No original atlas/JSON/XML/displacement files were changed.\n');

% Compare the original and manually refined VisuAlign registrations.
% Writes three images for each condition: histology, atlas, overlay.
% Does not modify the existing JSON, XML, or displacement CSV files.
clear; clc; close all;

%% Settings (edit these)
projectFolder = '/Users/minghanfan/Desktop/TongProject/TongProject_Code';
mouseID = 'CC1A';
sliceName = 'CC1A_163_A5_1';
lcID = 147;

mouseFolder = fullfile(projectFolder, 'rawData', mouseID);
infoFile = fullfile(mouseFolder, [mouseID '-info.xml']);
atlasFile = fullfile(projectFolder, 'atlas', 'annotationVolume.mat');

originalJson = fullfile(mouseFolder, [mouseID '-visualign.json']);
% Change this if you saved your manually refined JSON elsewhere.
manualJson = fullfile(mouseFolder, [mouseID '-visualign-manual.json']);

outputFolder = fullfile(projectFolder, 'LC_analysis', mouseID, ...
    'comparison', sliceName);
resizeFactor = 0.2;

%% Validate inputs
addpath(genpath(fullfile(projectFolder, 'brainAlignment')));
assert(isfile(infoFile), 'Info XML not found: %s', infoFile);
assert(isfile(atlasFile), 'Atlas MAT not found: %s', atlasFile);
assert(isfile(originalJson), 'Original JSON not found: %s', originalJson);
assert(isfile(manualJson), 'Manual JSON not found: %s', manualJson);

%% Load slice, original JSON record, and manual JSON record
sBase = Slice(sliceName, infoFile);
assert(~isempty(sBase.thumbnail), 'Missing thumbnail for %s.', sliceName);
recordOriginal = findSliceRecord(originalJson, sliceName);
recordManual = findSliceRecord(manualJson, sliceName);

% Set QuickNII planes explicitly using the corresponding JSON versions.
sOriginal = setRegistrationFromRecord(sBase, recordOriginal);
sManual = setRegistrationFromRecord(sBase, recordManual);

% The dataset's existing CSV displacement field belongs to the original
% registration. Keep it unchanged to reproduce the baseline overlays.
if isempty(sBase.dispField)
    warning('Original dispFields CSV files missing. Recomputing from original JSON markers.');
    sOriginal.dispField = calculateDisplacement(sOriginal, resizeFactor);
else
    sOriginal.dispField = sBase.dispField;
end

% The manual JSON has new marker positions. Recompute a displacement field
% from those markers: Slice does not do this automatically when loading JSON.
sManual.dispField = calculateDisplacement(sManual, resizeFactor);

%% Load and orient the atlas (same transformation as inspect_LC.m)
S = load(atlasFile, 'annotationVolume');
annotationVolume = permute(S.annotationVolume, [3 1 2]);
annotationVolume = flip(annotationVolume, 2);
annotationVolume = flip(annotationVolume, 3);

%% Apply each registration separately
planeOriginal = sOriginal.getObliqueSlice(annotationVolume);
atlasOriginal = imwarp(planeOriginal, sOriginal.dispField, 'nearest');

planeManual = sManual.getObliqueSlice(annotationVolume);
atlasManual = imwarp(planeManual, sManual.dispField, 'nearest');

% Defensive check: do not silently resize/misalign annotations.
imageSize = size(sBase.thumbnail, [1 2]);
assert(isequal(size(atlasOriginal), imageSize), ...
    'Original atlas image size differs from thumbnail size.');
assert(isequal(size(atlasManual), imageSize), ...
    'Manual atlas image size differs from thumbnail size.');

%% Save the six images
folderOriginal = fullfile(outputFolder, '01_Original');
folderManual = fullfile(outputFolder, '02_Manual_VisuAlign');
if ~isfolder(folderOriginal); mkdir(folderOriginal); end
if ~isfolder(folderManual); mkdir(folderManual); end

saveThreeImages(sBase.thumbnail, atlasOriginal, lcID, folderOriginal);
saveThreeImages(sBase.thumbnail, atlasManual, lcID, folderManual);

%% Summary
fprintf('\nDone: %s\n', sliceName);
fprintf('Original LC pixels: %d\n', nnz(atlasOriginal == lcID));
fprintf('Manual LC pixels:   %d\n', nnz(atlasManual == lcID));
fprintf('Output folder: %s\n', outputFolder);
fprintf('Original files were not modified.\n');

%% Local functions
function record = findSliceRecord(jsonFile, sliceName)
    % VisuAlign JSONs sometimes contain Windows-style paths.
    txt = fileread(jsonFile);
    J = jsondecode(strrep(txt, '\', '/'));
    assert(isfield(J, 'slices'), 'No "slices" field in %s', jsonFile);
    filenames = string({J.slices.filename});
    matches = find(contains(filenames, sliceName));
    assert(numel(matches) == 1, ...
        'Expected one record for %s in %s; found %d.', ...
        sliceName, jsonFile, numel(matches));
    record = J.slices(matches);
end

function s = setRegistrationFromRecord(s, record)
    a = double(record.anchoring(:)');
    assert(numel(a) == 9, 'Expected nine anchoring values.');
    s.o = a(1:3);
    s.u = a(4:6);
    s.v = a(7:9);
    s.markers = double(record.markers);
    assert(size(s.markers, 2) == 4, ...
        'Expected a markers array with four columns.');
end

function D = calculateDisplacement(s, resizeFactor)
    assert(~isempty(s.markers), ...
        'No VisuAlign markers in manually saved JSON for %s.', s.name);
    Dsmall = s.getDispField(resizeFactor);
    D = imresize(Dsmall, size(s.thumbnail, [1 2]), 'bilinear');
end

function saveThreeImages(histology, alignedAtlas, lcID, outputFolder)
    % 1: Original histology thumbnail, unchanged.
    rgb = histology;
    if size(rgb, 3) == 1
        rgb = repmat(rgb, [1 1 3]);
    end
    if ~isa(rgb, 'uint8')
        rgb = im2uint8(rgb);
    end
    imwrite(rgb, fullfile(outputFolder, 'histology.png'));

    % Anatomical boundaries from the label image (use same method as Slice.m).
    borders = Slice.annotationBorders(Slice.rescaleAnnotation(alignedAtlas));
    lcBoundary = bwperim(alignedAtlas == lcID, 8);
    lcFill = alignedAtlas == lcID;

    % 2: Atlas only. Dark background, light anatomical boundaries, red LC.
    atlasView = uint8(25 * ones([size(alignedAtlas,1),size(alignedAtlas,2),3]));
    atlasView = colorPixels(atlasView, lcFill, [110 30 30]);
    atlasView = colorPixels(atlasView, borders, [215 215 215]);
    atlasView = colorPixels(atlasView, lcBoundary, [255 55 55]);
    imwrite(atlasView, fullfile(outputFolder, 'atlas.png'));

    % 3: Overlay, white anatomical boundaries + red LC boundary.
    overlay = blendPixels(rgb, borders, [255 255 255], 0.50);
    overlay = blendPixels(overlay, lcBoundary, [255 30 30], 1.00);
    imwrite(overlay, fullfile(outputFolder, 'overlay.png'));
end

function img = colorPixels(img, mask, color)
    for c = 1:3
        ch = img(:,:,c);
        ch(mask) = uint8(color(c));
        img(:,:,c) = ch;
    end
end

function img = blendPixels(img, mask, color, alpha)
    img = double(img);
    for c = 1:3
        ch = img(:,:,c);
        ch(mask) = (1-alpha) * ch(mask) + alpha * color(c);
        img(:,:,c) = ch;
    end
    img = uint8(round(img));
end

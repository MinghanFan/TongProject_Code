% Find mouse slices whose reference atlas contains LC (ID 147),
% then save an LC overlay and an anatomical-boundary overlay.
clear; clc; close all;

%% Settings
mouseID = 'CC4A';
lcID = 147;
projectFolder = '/Users/minghanfan/Desktop/TongProject/TongProject_Code';

%% Paths
addpath(genpath(fullfile(projectFolder, 'brainAlignment')));
infoFile = fullfile(projectFolder, 'rawData', mouseID, [mouseID '-info.xml']);
atlasFile = fullfile(projectFolder, 'atlas', 'annotationVolume.mat');
outputFolder = fullfile(projectFolder, 'LC_analysis', mouseID);
imageFolder = fullfile(outputFolder, 'images');

if ~isfolder(imageFolder)
    mkdir(imageFolder);
end

%% Atlas orientation (retain the validated transformation)
S = load(atlasFile);
annotationVolume = permute(S.annotationVolume, [3 1 2]);
annotationVolume = flip(annotationVolume, 2);
annotationVolume = flip(annotationVolume, 3);

%% Find LC-containing reference sections and inspect their alignment
sliceArray = allSlicesFromXml(infoFile);
sliceNames = strings(0, 1);

fprintf('Searching %d slices for mouse %s...\n', numel(sliceArray), mouseID);

for i = 1:numel(sliceArray)
    s = sliceArray(i);

    % Skip slices without a thumbnail or atlas anchoring information.
    if isempty(s.thumbnail) || isempty(s.o) || ...
       isempty(s.u) || isempty(s.v)
        continue
    end

    % Selection is based ONLY on the unwarped reference atlas labels.
    obliqueSlice = s.getObliqueSlice(annotationVolume);
    if ~any(obliqueSlice(:) == lcID)
        continue
    end

    fprintf('\n%s: Reference atlas contains LC\n', s.name);
    sliceNames(end + 1, 1) = string(s.name);

    % VisuAlign is used for inspection images, not for slice selection.
    if isempty(s.dispField)
        fprintf('No VisuAlign displacement for %s. Skipping images.\n', s.name);
        continue
    end

    annot = imwarp(obliqueSlice, s.dispField, 'nearest');
    annot = imresize(annot, size(s.thumbnail, [1 2]), 'nearest');
    lcMask = (annot == lcID);

    % Image 1: LC boundary (red) on microscopy thumbnail.
    fig = figure('Visible', 'off');
    imshow(s.thumbnail);
    hold on;
    if any(lcMask(:))
        visboundaries(lcMask, 'Color', 'r');
    end
    title(['LC alignment: ' s.name], 'Interpreter', 'none');
    exportgraphics(fig, fullfile(imageFolder, [s.name '_LC_overlay.png']), ...
        'Resolution', 300);
    close(fig);

    % Image 2: All atlas boundaries (white) with LC boundary (red).
    rescaled = Slice.rescaleAnnotation(annot);
    borders = Slice.annotationBorders(rescaled);

    fig = figure('Visible', 'off');
    imshow(s.thumbnail);
    hold on;
    whiteImg = uint8(255 * ones(size(s.thumbnail, 1), size(s.thumbnail, 2)));
    overlay = imshow(whiteImg);
    overlay.AlphaData = double(borders) * 0.5;
    if any(lcMask(:))
        visboundaries(lcMask, 'Color', 'r');
    end
    title(['Atlas boundaries: ' s.name], 'Interpreter', 'none');
    exportgraphics(fig, fullfile(imageFolder, [s.name '_boundaries.png']), ...
        'Resolution', 300);
    close(fig);

    fprintf('Inspection images saved for %s\n', s.name);
end

%% Save CSV of all reference-atlas LC slices
results = table(sliceNames, 'VariableNames', {'SliceName'});
csvFile = fullfile(outputFolder, [mouseID '_LC_reference_slices.csv']);
writetable(results, csvFile);

fprintf('\nSearch complete! Mouse: %s\n', mouseID);
fprintf('Atlas sections containing LC: %d\n', height(results));
fprintf('CSV: %s\nImages: %s\n', csvFile, imageFolder);
disp(results);

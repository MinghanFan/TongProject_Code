% Inspect LC and atlas boundaries for manually selected mouse slices.
clear; clc; close all;

%% Settings
mouseID = 'CC1A';
sliceNames = {
    'CC1A_070_A4_1'
};
lcID = 147;
projectFolder = '/Users/minghanfan/Desktop/TongProject/TongProject_Code';

%% Paths
addpath(genpath(fullfile(projectFolder, 'brainAlignment')));
infoFile = fullfile(projectFolder, 'rawData', mouseID, [mouseID '-info.xml']);
atlasFile = fullfile(projectFolder, 'atlas', 'annotationVolume.mat');
outputFolder = fullfile(projectFolder, 'LC_analysis');

if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

%% Atlas orientation (retain the validated transformation)
S = load(atlasFile);
annotationVolume = permute(S.annotationVolume, [3 1 2]);
annotationVolume = flip(annotationVolume, 2);
annotationVolume = flip(annotationVolume, 3);

%% Inspect selected slices
for i = 1:numel(sliceNames)
    s = Slice(sliceNames{i}, infoFile);
    fprintf('\nProcessing %s\n', s.name);

    % Extract reference atlas section.
    obliqueSlice = s.getObliqueSlice(annotationVolume);
    coverage = 100 * nnz(obliqueSlice) / numel(obliqueSlice);
    fprintf('Atlas coverage: %.2f%%\n', coverage);
    fprintf('LC pixels before warp: %d\n', nnz(obliqueSlice == lcID));

    % Apply saved VisuAlign displacement to the annotation labels.
    assert(~isempty(s.dispField), 'Displacement field not loaded');
    annot = imwarp(obliqueSlice, s.dispField, 'nearest');
    annot = imresize(annot, size(s.thumbnail, [1 2]), 'nearest');
    lcMask = (annot == lcID);
    fprintf('LC pixels after warp: %d\n', nnz(lcMask));

    % Image 1: LC boundary (red) on microscopy thumbnail.
    fig = figure('Visible', 'off');
    imshow(s.thumbnail);
    hold on;
    if any(lcMask(:))
        visboundaries(lcMask, 'Color', 'r');
    end
    title(['LC alignment: ' s.name], 'Interpreter', 'none');
    exportgraphics(fig, fullfile(outputFolder, [s.name '_LC_overlay.png']));
    close(fig);

    % Image 2: All atlas boundaries (white) with LC boundary (red).
    rescaled = Slice.rescaleAnnotation(annot);
    borders = Slice.annotationBorders(rescaled);

    fig = figure('Visible', 'off');
    imshow(s.thumbnail);
    hold on;
    whiteImg = uint8(255 * ones(size(s.thumbnail)));
    overlay = imshow(whiteImg);
    overlay.AlphaData = double(borders) * 0.5;
    if any(lcMask(:))
        visboundaries(lcMask, 'Color', 'r');
    end
    title(['Atlas boundaries: ' s.name], 'Interpreter', 'none');
    exportgraphics(fig, fullfile(outputFolder, [s.name '_boundaries1.png']));
    close(fig);
end

fprintf('Finished! Images saved to: %s\n', outputFolder);

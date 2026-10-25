clear; clc; close all;

% Set up paths
projectFolder = '/Users/minghanfan/Desktop/TongProject/TongProject_Code';
addpath(genpath(fullfile(projectFolder, 'brainAlignment')));

infoFile = fullfile(projectFolder, 'rawData', 'CC1A', 'CC1A-info.xml');
atlasFile = fullfile(projectFolder, 'atlas', 'annotationVolume.mat');
outputFolder = fullfile(projectFolder, 'LC_analysis');

% Load and reorder Allen CCFv3 annotation
S = load(atlasFile);
annotationVolume = permute(S.annotationVolume, [3 1 2]);
annotationVolume = flip(annotationVolume, 2);
annotationVolume = flip(annotationVolume, 3);

% LC region ID
lcID = 147;

% Slices containing LC in the atlas
sliceNames = {
    % 'CC1A_023_B2_1'
    % 'CC1A_070_A4_1'
    % 'CC1A_091_A1_1'
    % 'CC1A_154_A2_1'
    % 'CC1A_157_A3_1'
    % 'CC1A_160_A4_1'
    % 'CC1A_163_A5_1'
    % 'CC1A_166_A6_1'
    % 'CC1A_070_A4_1'
    % 'CC4A_007_C3_2'
    % 'CC4A_040_C4_1'
    % 'CC2B_001_A1_1'
    % 'CC2B_082_A1_1'
    % 'AL1A_136_A1_1'
    % 'AL1A_157_B2_1'
    % 'CC3A_121_A5_2'
};

% Process each slice
for i = 1:numel(sliceNames)

    s = Slice(sliceNames{i}, infoFile);
    fprintf('\nProcessing %s\n', s.name);

    % Extract atlas section using QuickNII alignment
    obliqueSlice = s.getObliqueSlice(annotationVolume);

    coverage = 100 * nnz(obliqueSlice) / numel(obliqueSlice);
    fprintf('Atlas coverage: %.2f%%\n', coverage);
    fprintf('LC pixels before warp: %d\n', nnz(obliqueSlice == lcID));

    % Apply VisuAlign displacement
    assert(~isempty(s.dispField), 'Displacement field not loaded');

    annot = imwarp(obliqueSlice, s.dispField, 'nearest');
    annot = imresize(annot, size(s.thumbnail, [1 2]), 'nearest');

    % Identify LC
    lcMask = annot == lcID;
    fprintf('LC pixels after warp: %d\n', nnz(lcMask));

    % Figure 1: LC boundary
    fig = figure('Visible', 'off');
    imshow(s.thumbnail);
    hold on;

    if any(lcMask(:))
        visboundaries(lcMask, 'Color', 'r');
    end

    title(['LC alignment: ' s.name]);
    exportgraphics(fig, fullfile(outputFolder, ...
        [s.name '_LC_overlay.png']));
    close(fig);

    % Figure 2: All anatomical boundaries and LC
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

    title(['Atlas boundaries: ' s.name]);
    exportgraphics(fig, fullfile(outputFolder, ...
        [s.name '_boundaries.png']));
    close(fig);

end

disp('Finished! All images saved to LC_analysis.');

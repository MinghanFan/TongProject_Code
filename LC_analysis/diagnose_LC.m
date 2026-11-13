% Diagnose how many QuickNII-aligned atlas planes intersect LC (ID 147).
% LC is tested ONLY in the extracted atlas, not in the microscopy image.
% VisuAlign is used ONLY for visual inspection, never to select slices.
clear; clc; close all;

%% Settings
mouseID = 'CC1B';   % Change to CC1A, CC3A, CC4B, etc.
lcID = 147;
voxelSizeUm = 25;
nearRadiusVox = 4;  % 4 x 25 = 100 um spatial proximity (not LC itself)
saveImages = true;
projectFolder = '/Users/minghanfan/Desktop/TongProject/TongProject_Code';

%% Paths
addpath(genpath(fullfile(projectFolder, 'brainAlignment')));
infoFile = fullfile(projectFolder, 'rawData', mouseID, [mouseID '-info.xml']);
atlasFile = fullfile(projectFolder, 'atlas', 'annotationVolume.mat');
outFolder = fullfile(projectFolder, 'LC_analysis', mouseID, 'LC_diagnostic');
imgFolder = fullfile(outFolder, 'images');
if ~isfolder(imgFolder), mkdir(imgFolder); end

%% Load atlas using your validated orientation (do not change flips)
S = load(atlasFile);
annotationVolume = permute(S.annotationVolume, [3 1 2]);
annotationVolume = flip(annotationVolume, 2);
annotationVolume = flip(annotationVolume, 3);
lc3D = (annotationVolume == lcID);

% From the QuickNII anchoring-vector orientation used in this project:
% dimension 1 = left/right (ML), dimension 2 = AP, dimension 3 = DV.
% Check these directions against your reference atlas if it differs.
midML = floor(size(lc3D, 1)/2);
side1 = reshape(sum(sum(lc3D(1:midML,:,:), 1), 3), [], 1);
side2 = reshape(sum(sum(lc3D(midML+1:end,:,:), 1), 3), [], 1);
apIndices = find(side1 + side2 > 0);
if isempty(apIndices), error('LC (ID 147) not found in atlas.'); end

fprintf('\n=== ORIGINAL 3D ATLAS ===\n');
fprintf('Atlas dimensions: %s\n', mat2str(size(annotationVolume)));
fprintf('LC voxels total: %d\n', nnz(lc3D));
fprintf('ML half totals: %d vs %d voxels\n', sum(side1), sum(side2));
fprintf('LC AP index range: %d to %d\n', min(apIndices), max(apIndices));
fprintf('AP planes with LC: %d\n', numel(apIndices));
fprintf('AP extent (range): %.0f um\n', ...
    (max(apIndices)-min(apIndices)+1)*voxelSizeUm);
fprintf('Both sides: %d AP planes\n', nnz(side1>0 & side2>0));
fprintf('Only side 1: %d AP planes\n', nnz(side1>0 & side2==0));
fprintf('Only side 2: %d AP planes\n', nnz(side1==0 & side2>0));

%% Create a nearby-LC region for DIAGNOSTICS ONLY
% Includes locations within a 3D radius of 100 um of atlas LC.
% It is NOT an alternative LC annotation.
[r1,r2,r3] = ndgrid(-nearRadiusVox:nearRadiusVox, ...
                        -nearRadiusVox:nearRadiusVox, ...
                        -nearRadiusVox:nearRadiusVox);
ball = (r1.^2 + r2.^2 + r3.^2) <= nearRadiusVox^2;
near3D = uint8(imdilate(lc3D, ball));

%% Read all microscopy sections and test their QuickNII atlas planes
sliceArray = allSlicesFromXml(infoFile);
n = numel(sliceArray);
sliceName = strings(n,1);
apCenter = nan(n,1);
lcPixels = nan(n,1);
nearPixels = nan(n,1);
status = strings(n,1);

fprintf('\n=== SEARCHING %d SECTIONS: %s ===\n', n, mouseID);
for i = 1:n
    s = sliceArray(i);
    sliceName(i) = string(s.name);

    if isempty(s.thumbnail) || isempty(s.o) || isempty(s.u) || isempty(s.v)
        status(i) = "Missing image or QuickNII anchors";
        continue
    end

    % Approximate AP coordinate at the center of the QuickNII plane.
    apCenter(i) = s.o(2) + 0.5*s.u(2) + 0.5*s.v(2);

    try
        % Exact LC: original selection rule, with no VisuAlign involved.
        plane = s.getObliqueSlice(annotationVolume);
        lcPixels(i) = nnz(plane == lcID);

        % Near LC: same QuickNII plane, but intersecting the 100 um envelope.
        nearPlane = s.getObliqueSlice(near3D);
        nearPixels(i) = nnz(nearPlane > 0);
    catch ME
        status(i) = "Atlas extraction error";
        fprintf('ERROR: %s: %s\n', s.name, ME.message);
        continue
    end

    if lcPixels(i) > 0
        status(i) = "Exact LC";
    elseif nearPixels(i) > 0
        status(i) = "Near LC (not exact)";
    else
        status(i) = "No LC";
    end

    fprintf('%-25s AP center %7.1f | exact %6d | near %6d | %s\n', ...
        s.name, apCenter(i), lcPixels(i), nearPixels(i), status(i));

    % Inspection images for exact hits and near misses.
    if ~saveImages || nearPixels(i) == 0, continue; end
    targetSize = size(s.thumbnail, [1 2]);
    rawLC = imresize(plane == lcID, targetSize, 'nearest');
    rawNear = imresize(nearPlane > 0, targetSize, 'nearest');

    haveDisp = ~isempty(s.dispField);
    if haveDisp
        warped = imwarp(plane, s.dispField, 'nearest');
        warpedNear = imwarp(nearPlane > 0, s.dispField, 'nearest');
        postLC = imresize(warped == lcID, targetSize, 'nearest');
        postNear = imresize(warpedNear > 0, targetSize, 'nearest');
    end

    fig = figure('Visible','off', 'Position',[100 100 1300 600]);
    tiledlayout(1,2);
    nexttile; imshow(s.thumbnail); hold on;
    if any(rawNear(:)), visboundaries(rawNear,'Color','y'); end
    if any(rawLC(:)), visboundaries(rawLC,'Color','r'); end
    title('QuickNII only', 'Interpreter','none');

    nexttile; imshow(s.thumbnail); hold on;
    if haveDisp
        if any(postNear(:)), visboundaries(postNear,'Color','y'); end
        if any(postLC(:)), visboundaries(postLC,'Color','r'); end
        title('After VisuAlign (inspection only)');
    else
        title('No VisuAlign field available');
    end
    sgtitle(sprintf('%s | red=actual LC, yellow=within 100 um', s.name), ...
        'Interpreter','none');
    exportgraphics(fig, fullfile(imgFolder, [char(sliceName(i)) '_LC_diagnostic.png']), ...
        'Resolution',150);
    close(fig);
end

%% Save results for ALL sections, including those without LC
results = table(sliceName, apCenter, lcPixels, nearPixels, status, ...
    'VariableNames', {'SliceName','QuickNII_AP_Center_Index', ...
                      'Exact_LC_Pixels','Within_100um_Pixels','Status'});
csvFile = fullfile(outFolder, [mouseID '_LC_diagnostic.csv']);
writetable(results, csvFile);

%% Plot the corrected coronal AP distribution and section locations
fig = figure('Visible','off','Position',[100 100 1200 550]);
plot(side1, 'b-', 'LineWidth',1.5); hold on;
plot(side2, 'r-', 'LineWidth',1.5);
yTop = max([side1;side2]);
yMarker = yTop*1.08;
exact = lcPixels > 0;
nearOnly = lcPixels == 0 & nearPixels > 0;
other = lcPixels == 0 & nearPixels == 0;
plot(apCenter(exact), yMarker*ones(nnz(exact),1), 'kv', 'MarkerFaceColor','k');
plot(apCenter(nearOnly), yMarker*ones(nnz(nearOnly),1), 'mv', 'MarkerFaceColor','m');
plot(apCenter(other), zeros(nnz(other),1), 'k.', 'MarkerSize',7);
xlabel('Atlas AP index (dim 2)'); ylabel('LC voxels per AP plane');
title([mouseID ' | CCF LC distribution and QuickNII section positions']);
legend('Hemisphere 1','Hemisphere 2','Exact LC sections', ...
       'Near misses (100 um)','Other sections', 'Location','best');
ylim([0 max(1,yTop*1.22)]); grid on;
exportgraphics(fig, fullfile(outFolder, [mouseID '_LC_AP_distribution.png']), ...
    'Resolution',180);
close(fig);

fprintf('\n=== SUMMARY: %s ===\n', mouseID);
fprintf('Total sections: %d\n', n);
fprintf('Valid QuickNII sections: %d\n', nnz(isfinite(lcPixels)));
fprintf('Exact LC matches: %d\n', nnz(exact));
fprintf('Near misses (within 100 um): %d\n', nnz(nearOnly));
fprintf('Other sections: %d\n', nnz(other));
fprintf('CSV: %s\n', csvFile);
fprintf('Images: %s\n', imgFolder);
disp(results(exact | nearOnly, :));

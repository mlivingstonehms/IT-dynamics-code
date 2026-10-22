%% calculate_spatial_frequency_manyOO_cpd.m
%
% Calculate a spatial-frequency metric for every manyOO image,
% in EXACTLY the same image-name order as manyOOimagename.mat.
%
% Main output:
%   spatialFreqCpd = power-weighted radial spectral centroid, cycles/degree
%
% Assumption:
%   every stimulus spans 5 degrees horizontally.
%
% Improvements over the older FFT code:
%   - removes mean luminance (DC)
%   - applies a 2-D Hann window to reduce image-edge artifacts
%   - explicitly removes the zero-frequency bin
%   - expresses frequency directly in cycles/degree
%
% IMPORTANT:
%   Analyze the ENTIRE displayed image, including white background.
%   Do not apply the non-white mask used for colorfulness, because the
%   edge of such a mask would itself add artificial high-frequency power.

close all
clearvars
clc

%% PATHS

imageFolder = ...
    '/n/data2/hms/neurobio/livingstone/Stimuli/manyOO/';

metricFolder = ...
    '/n/data2/hms/neurobio/livingstone/marge/margemonkeys/chatmetricsEdges/';

nameFile = ...
    fullfile(metricFolder,'manyOOimagename.mat');

saveFile = ...
    fullfile(metricFolder,'spatialFrequencyCpd_manyOO.mat' ...
    );
%% SETTINGS

imageWidthDeg = 5;
checkpointEvery = 250;

%% LOAD MASTER IMAGE ORDER

S = load(nameFile);

if ~isfield(S,'imageName')
    error('Expected variable "imageName" in %s',nameFile);
end

imageName = string(S.imageName(:));
nImages = numel(imageName);

fprintf('Master manyOO image list contains %d images.\n',nImages);

%% BUILD RECURSIVE FILE LOOKUP

extensions = { ...
    '*.jpg','*.JPG', ...
    '*.jpeg','*.JPEG', ...
    '*.png','*.PNG', ...
    '*.tif','*.TIF', ...
    '*.tiff','*.TIFF'};

allFiles = [];

for e = 1:numel(extensions)
    theseFiles = dir(fullfile(imageFolder,'**',extensions{e}));
    allFiles = [allFiles; theseFiles]; %#ok<AGROW>
end

fprintf('Found %d image files under imageFolder.\n',numel(allFiles));

diskNames = strings(numel(allFiles),1);
diskPaths = strings(numel(allFiles),1);

for k = 1:numel(allFiles)
    diskNames(k) = lower(string(allFiles(k).name));
    diskPaths(k) = string(fullfile(allFiles(k).folder,allFiles(k).name));
end

%% CHECK FOR DUPLICATE FILENAMES

[uNames,~,ic] = unique(diskNames);
counts = accumarray(ic,1);
duplicateNames = uNames(counts > 1);

if ~isempty(duplicateNames)
    fprintf('\nWARNING: %d filenames occur more than once under manyOO.\n', ...
        numel(duplicateNames));
    fprintf('The first matching file will be used for duplicate basenames.\n\n');
end

%% PREALLOCATE

spatialFreqCpd = nan(nImages,1);
spectralSpreadCpd = nan(nImages,1);
totalSpectralPower = nan(nImages,1);
imageFound = false(nImages,1);

%% LOOP THROUGH IMAGES IN MASTER imageName ORDER

for im = 1:nImages

    targetName = lower(imageName(im));
    idx = find(diskNames == targetName,1,'first');

    if isempty(idx)
        warning('Could not find image %d/%d: %s', ...
            im,nImages,imageName(im));
        continue
    end

    imageFound(im) = true;

    %% Read image

    [I,map] = imread(diskPaths(idx));

    if ~isempty(map)
        I = ind2rgb(I,map);
    end

    % Convert RGB to grayscale
    if ndims(I) == 3
        I = I(:,:,1:min(3,size(I,3)));

        if size(I,3) == 3
            I = rgb2gray(I);
        else
            I = I(:,:,1);
        end
    end

    I = double(I);

    [M,N] = size(I);

    if M < 2 || N < 2
        warning('Image too small for FFT analysis: %s',imageName(im));
        continue
    end

    %% Remove mean luminance / DC

    I(~isfinite(I)) = 0;
    I = I - mean(I(:));

    %% Apply 2-D Hann window
    % Written explicitly so hann() is not required.

    wy = 0.5 - 0.5*cos(2*pi*(0:M-1)'/(M-1));
    wx = 0.5 - 0.5*cos(2*pi*(0:N-1)/(N-1));

    window2D = wy * wx;

    Iw = I .* window2D;

    %% Fourier power spectrum

    F = fftshift(fft2(Iw));
    P = abs(F).^2;

    %% Frequency coordinates in cycles/degree
    %
    % Width = 5 deg.
    % For a non-square image, infer vertical visual angle from aspect ratio,
    % assuming square pixels.

    imageHeightDeg = imageWidthDeg * (M/N);

    fx = (-floor(N/2):ceil(N/2)-1) / imageWidthDeg;
    fy = (-floor(M/2):ceil(M/2)-1) / imageHeightDeg;

    [FX,FY] = meshgrid(fx,fy);

    radialFreqCpd = sqrt(FX.^2 + FY.^2);

    %% Remove the zero-frequency bin explicitly

    P(radialFreqCpd == 0) = 0;

    totalPower = sum(P(:));

    if totalPower <= 0
        warning('No non-DC spectral power for image: %s',imageName(im));
        continue
    end

    %% MAIN METRIC: power-weighted radial spectral centroid

    sf = sum(radialFreqCpd(:).*P(:)) / totalPower;

    spatialFreqCpd(im) = sf;
    totalSpectralPower(im) = totalPower;

    %% Optional diagnostic: spectral spread

    spectralSpreadCpd(im) = sqrt( ...
        sum(((radialFreqCpd(:)-sf).^2).*P(:)) / totalPower);

    %% Progress

    if mod(im,100)==0 || im==1 || im==nImages
        fprintf('%5d/%5d  %-35s  SF = %.3f cycles/deg\n', ...
            im,nImages,imageName(im),spatialFreqCpd(im));
    end

    %% Checkpoint

    if mod(im,checkpointEvery)==0

        save(saveFile, ...
            'spatialFreqCpd', ...
            'spectralSpreadCpd', ...
            'totalSpectralPower', ...
            'imageFound', ...
            'imageName', ...
            'imageWidthDeg', ...
            '-v7.3');

        fprintf('Checkpoint saved at image %d.\n',im);
    end
end

%% RANK NORMALIZE

spatialFreqRank = nan(size(spatialFreqCpd));

valid = isfinite(spatialFreqCpd);

r = tiedrank(spatialFreqCpd(valid));

spatialFreqRank(valid) = ...
    (r - 0.5) / sum(valid);

%% FINAL SAVE

save(saveFile, ...
    'spatialFreqCpd', ...
    'spatialFreqRank', ...
    'spectralSpreadCpd', ...
    'totalSpectralPower', ...
    'imageFound', ...
    'imageName', ...
    'imageWidthDeg', ...
    '-v7.3');

fprintf('\nFinished.\n');
fprintf('Images successfully found: %d / %d\n',sum(imageFound),nImages);
fprintf('Valid spatial-frequency values: %d / %d\n',sum(valid),nImages);
fprintf('Saved:\n%s\n',saveFile);

%% QUICK SANITY-CHECK FIGURES

figure('Color','w');
histogram(spatialFreqCpd(valid),100);

xlabel('Spatial-frequency centroid (cycles/degree)');
ylabel('Number of images');
title('manyOO spatial-frequency distribution');

set(gca,'TickDir','out','LineWidth',1.5,'FontSize',12);
box off

figure('Color','w');
histogram(spatialFreqRank(valid),50);

xlabel('Rank-normalized spatial frequency');
ylabel('Number of images');
title('Rank-normalized manyOO spatial frequency');

set(gca,'TickDir','out','LineWidth',1.5,'FontSize',12);
box off

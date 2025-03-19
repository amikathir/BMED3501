%% part3_analysis.m
% Updated Part 3: Spatial/Frequency Domain Analysis
% ------------------------------------------------
% Demonstrates:
%  1) Central Slice Theorem
%  2) Fourier coverage diagrams for multiple angular sampling sets
%  3) Aliasing artifact progression (optional sections)

clear; clc; close all;

%% 1. Generate a phantom and its 2D Fourier transform
N = 256; % Phantom/image size
phantomImg = phantom(N);

% 2D Fourier transform of the phantom
F2_phantom = fftshift(fft2(phantomImg));

% Frequency axes for visualization
freqAxis = linspace(-0.5, 0.5, N);
[Xf, Yf] = meshgrid(freqAxis, freqAxis);

% Plot original phantom and its 2D Fourier magnitude
figure('Name','Phantom and 2D Fourier Magnitude');
subplot(1,2,1);
imshow(phantomImg, []);
title('Original Phantom');
subplot(1,2,2);
imshow(log(abs(F2_phantom)+1), []);
title('2D FFT of Phantom (log scale)');
drawnow;


%% 2. Central Slice Theorem Demonstration
theta0 = 45; % pick one angle
[sinogramSingle, xp] = radon(phantomImg, theta0);

% The sinogram has just one column if we pass a single angle
projection = sinogramSingle(:,1);
nDet = length(projection);

% 1D FFT of that projection
projFFT = fftshift(fft(projection));
freqProj = linspace(-0.5, 0.5, nDet);

% Plot the radial slice in the 2D FFT domain
figure('Name','Central Slice Theorem');
imshow(log(abs(F2_phantom)+1), []);
title(sprintf('Central Slice at \\theta = %d^{\\circ}', theta0));
hold on;

% Draw a line corresponding to the radial slice
r = linspace(-0.5, 0.5, nDet);  % radial freq axis
xf_line = r * cosd(theta0);
yf_line = r * sind(theta0);

% Convert freq to pixel coords
xf_px = (xf_line + 0.5) * N; 
yf_px = (yf_line + 0.5) * N;
plot(xf_px, yf_px, 'r', 'LineWidth',1.5);

figure('Name','1D FFT of Projection');
plot(freqProj, abs(projFFT), 'LineWidth',1.5);
title('Magnitude of 1D FFT of Projection');
xlabel('Normalized Frequency'); ylabel('|FFT|');
grid on;
drawnow;


%% 3. Fourier Coverage for Different Angular Sampling
% We expand anglesSet to include multiple increments, e.g. 1°, 2°, 5°, 10°.
anglesSet = {
    0:1:179,  ...
    0:2:179,  ...
    0:5:179,  ...
    0:10:179
};

figure('Name','Fourier Coverage for Multiple Theta Samplings');
colormap('jet');

numPlots = length(anglesSet);
for i = 1:numPlots
    angles = anglesSet{i};
    % Generate sinogram for each set of angles
    sino = radon(phantomImg, angles);
    [numDet, numAng] = size(sino);
    
    % Plot log(2D FFT) in the background
    subplot(1, numPlots, i);
    imagesc(log(abs(F2_phantom)+1)); axis image;
    hold on; title(sprintf('\\Delta\\theta = %d^{\\circ}', angles(2)-angles(1)));
    
    % Plot coverage lines for each angle
    for k = 1:numAng
        theta_k = angles(k);
        freq_s  = linspace(-0.5, 0.5, numDet);
        xf = freq_s .* cosd(theta_k);
        yf = freq_s .* sind(theta_k);
        
        % Convert freq coords to pixel coords
        xf_px = (xf + 0.5) * N;
        yf_px = (yf + 0.5) * N;
        plot(xf_px, yf_px, 'w');  % white lines over log-FFT
    end
    
    axis off;
end
sgtitle('Comparison of Fourier Coverage for Different Theta Increments');
drawnow;

%% Aliasing Artifact Progression with Zoomed-In Section
% (1°, 2°, 5°, 10°)

% 1. Define the sets of angles
anglesSet = {
    0:1:179,  ...  % 1-degree steps
    0:2:179,  ...  % 2-degree steps
    0:5:179,  ...  % 5-degree steps
    0:10:179  ...  % 10-degree steps
};

% 2. Generate a Shepp-Logan phantom (if not already done)
N = 256;
phantomImg = phantom(N);

% 3. Prepare figure
figure('Name','Aliasing Artifacts with Zoomed View');
numSets = length(anglesSet);

% Define a central ROI for zooming
roiSize = 50;                    % side length of the zoomed area
centerCoord = N/2;               % approximate center of the image
roiStart = round(centerCoord - roiSize/2) + 1;
roiEnd   = roiStart + roiSize - 1;

for i = 1:numSets
    angles = anglesSet{i};
    samplingStep = angles(2) - angles(1);
    
    % Create the sinogram
    sino = radon(phantomImg, angles);
    
    % Reconstruct using built-in FBP with Ram-Lak filter
    recon = iradon(sino, angles, 'linear','Ram-Lak', 1.0, N);

    % Subplot for full reconstruction
    subplot(2, numSets, i);
    imshow(recon, []);
    title(sprintf('Full View (%d\\circ steps)', samplingStep));

    % Subplot for zoomed region
    subplot(2, numSets, i + numSets);
    zoomedRegion = recon(roiStart:roiEnd, roiStart:roiEnd);
    imshow(zoomedRegion, []);
    title(sprintf('Zoomed (%d\\circ)', samplingStep));
end

sgtitle('Aliasing Artifacts: Angular Undersampling (Full vs. Zoomed Views)');
drawnow;

%% Aliasing Artifact Progression: Each Angular Range in Its Own Figure
% (1°, 2°, 5°, 10°)

% 1. Define the sets of angles
anglesSet = {
    0:1:179,  ...  % 1-degree steps
    0:2:179,  ...  % 2-degree steps
    0:5:179,  ...  % 5-degree steps
    0:10:179  ...  % 10-degree steps
};

% 2. Generate a Shepp-Logan phantom (if not already done)
N = 256;
phantomImg = phantom(N);

% Define a central ROI for zooming
roiSize = 50;                   
centerCoord = N/2;             
roiStart = round(centerCoord - roiSize/2) + 1;
roiEnd   = roiStart + roiSize - 1;

% 3. Loop over each angular set, reconstruct, and open two figures per set
for i = 1:length(anglesSet)
    angles = anglesSet{i};
    samplingStep = angles(2) - angles(1);
    
    % Create the sinogram
    sino = radon(phantomImg, angles);
    
    % Reconstruct using built-in FBP with Ram-Lak filter
    recon = iradon(sino, angles, 'linear','Ram-Lak', 1.0, N);
    
    % Zoomed region
    zoomedRegion = recon(roiStart:roiEnd, roiStart:roiEnd);
    
    %% Show the full reconstruction in its own figure
    figure('Name', sprintf('Full Reconstruction: %d-degree steps', samplingStep));
    imshow(recon, []);
    title(sprintf('Full View (%d\\circ steps)', samplingStep));
    drawnow;
    
    %% Show the zoomed region in a separate figure
    figure('Name', sprintf('Zoomed Reconstruction: %d-degree steps', samplingStep));
    imshow(zoomedRegion, []);
    title(sprintf('Zoomed View (%d\\circ steps)', samplingStep));
    drawnow;
end


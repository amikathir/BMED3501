clear; clc; close all;

phantomSize = 256;       % Image size
numAngles = 180;         % Number of projection angles
theta = linspace(0, 179, numAngles);

phantomImg = phantom(phantomSize);

[sinogram, xp] = radon(phantomImg, theta);

figure;
subplot(1,2,1);
imshow(phantomImg, []);
title('Original Phantom');

subplot(1,2,2);
imagesc(theta, xp, sinogram);
colormap gray; colorbar;
title('Sinogram');
xlabel('Angle (degrees)'); ylabel('Detector Position');

reconImg = FFBP(sinogram, theta);

figure;
imshow(reconImg, []);
title('Reconstructed Image (FFBP)');

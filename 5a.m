clear; clc; close all;

disp('--- 5.1 High-Frequency Noise Amplification ---');

N = 128; % keep it modest for faster iterative reconstruction
phantomImg = phantom(N);
theta = 0:2:178; % step of 2 degrees
[cleanSino, xp] = radon(phantomImg, theta);

rng(0);
noiseLevel = 0.02;
noisySino = cleanSino + noiseLevel*max(cleanSino(:))*randn(size(cleanSino));

recon_FFBP = iradon(noisySino, theta, 'linear','Ram-Lak', 1.0, N);

numIters = 20; % tweak as desired
recon_Iter = iterativeSIRT(noisySino, theta, N, numIters);

figure('Name','5.1: High-Frequency Noise Amplification');
subplot(1,2,1);
imshow(recon_FFBP, []);
title('FFBP (Ram-Lak)');
subplot(1,2,2);
imshow(recon_Iter, []);
title(sprintf('Simple Iterative (%d iters)', numIters));
sgtitle('Comparison: FFBP vs. Iterative Reconstruction');

disp('--- 5.2 Truncation Artifacts ---');

fullSino = radon(phantomImg, theta);
numDetFull = size(fullSino,1);

keepFraction = 0.6;
numDetTrunc = round(keepFraction * numDetFull);
midIdx = round(numDetFull/2);
startIdx = midIdx - round(numDetTrunc/2);
endIdx   = startIdx + numDetTrunc - 1;

truncSino = fullSino(startIdx:endIdx, :);

recon_Full = iradon(fullSino, theta, 'linear','Ram-Lak', 1.0, N);
recon_Trunc = iradon(truncSino, theta, 'linear','Ram-Lak', 1.0, N);

figure('Name','5.2: Truncation Artifacts');
subplot(1,2,1);
imshow(recon_Full, []);
title('Full Sinogram Reconstruction');
subplot(1,2,2);
imshow(recon_Trunc, []);
title('Truncated Sinogram Reconstruction');
sgtitle('Effect of Truncation: Notice Boundary Artifacts');

disp('--- 5.3 Gibbs Phenomenon ---');

stepImg = zeros(N);
stepImg(:, 1:N/2) = 1.0; % left half is 1, right half is 0

stepSino = radon(stepImg, theta);
recon_Step = iradon(stepSino, theta, 'linear','Ram-Lak', 1.0, N);

midRow = round(N/2);
profile_Recon = recon_Step(midRow, :);

figure('Name','5.3: Gibbs Phenomenon');
subplot(1,2,1);
imshow(stepImg, []);
title('Step Phantom');
subplot(1,2,2);
plot(profile_Recon, 'LineWidth',1.5);
xlabel('Pixel');
ylabel('Reconstructed Intensity');
title('Profile at Center Row (Gibbs Overshoot)');
sgtitle('Gibbs Phenomenon at Sharp Boundary');

disp('--- 5.4 Ring Artifacts ---');

phantomImg2 = phantom(N);
ringSino = radon(phantomImg2, theta);

dIdx = 60; % pick some row in the sinogram
ringOffset = 0.3 * max(ringSino(:));
ringSino(dIdx, :) = ringSino(dIdx, :) + ringOffset; % add offset

recon_Ring = iradon(ringSino, theta, 'linear','Ram-Lak', 1.0, N);

figure('Name','5.4: Ring Artifacts');
subplot(1,2,1);
imagesc(ringSino); colormap gray; colorbar;
title('Sinogram with "Bad" Detector Row');
xlabel('Angle Index'); ylabel('Detector Row');
subplot(1,2,2);
imshow(recon_Ring, []);
title('Reconstruction with Ring Artifact');
sgtitle('Ring Artifact from Detector Miscalibration');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------- Helper Functions --------------- %

function xRecon = iterativeSIRT(sinogram, theta, N, numIters)

xRecon = zeros(N, N);
alpha = 0.5;          
for i = 1:numIters
    estSino = radon(xRecon, theta);
    
    R = sinogram - estSino;
    
    corrImg = iradon(R, theta, 'linear','none', 1.0, N);
    
    xRecon = xRecon + alpha * corrImg;
end
end

clear; clc; close all;

N = 256;
phantomImg = phantom(N);
theta = 0:179;  % 180 projections
[cleanSino, xp] = radon(phantomImg, theta);

recon_clean = iradon(cleanSino, theta, 'linear','Ram-Lak', 1.0, N);

R_with_rings = cleanSino;

ringRows = round([N/4, N/2, 3*N/4]); % pick 3 rows in the mid/center region
numAngles = size(cleanSino,2);

ringAmplitude = 0.3 * max(cleanSino(:));  % Adjust as desired
ringPattern   = ringAmplitude * cos(2*pi*(1:numAngles)/20);

for row = ringRows
    R_with_rings(row, :) = R_with_rings(row, :) + ringPattern;
end

[numDet, ~] = size(cleanSino);
N_pad = 2^nextpow2(2*numDet);
halfFreq = (0:(N_pad/2 - 1))/N_pad;  % 0..0.5
ramp = [halfFreq, fliplr(halfFreq)];

freqSym = linspace(-0.5, 0.5, N_pad);
hannWin = hannWindow(freqSym, halfFreq(end));
hannFilter = ramp .* hannWin;

recon_with_rings_no_filter = iradon(R_with_rings, theta, 'linear','none', 1.0, N);

filteredSino = applyFFBPFilter(R_with_rings, hannFilter, N_pad);
recon_with_rings_filtered = iradon(filteredSino, theta, 'linear','none', 1.0, N);

figure('Name','Ring Artifact Demonstration');

subplot(1,3,1);
imshow(recon_clean, []);
title('Reconstruction without Ring Artifacts');

subplot(1,3,2);
imshow(recon_with_rings_no_filter, []);
title('Ring Artifact (No Filter)');

subplot(1,3,3);
imshow(recon_with_rings_filtered, []);
title('Ring Artifact (Hann Filter)');

drawnow;

%% ------------------ Helper Functions ------------------ %%
function y = hannWindow(freq, fcut)
    y = zeros(size(freq));
    mask = (abs(freq) <= abs(fcut));
    x = freq(mask)/fcut;
    y(mask) = 0.5*(1 + cos(pi*x));
end

function filteredSino = applyFFBPFilter(sinogram, filterCurve, N_pad)
    [numDet, numAng] = size(sinogram);
    filteredSino = zeros(size(sinogram));
    for angIdx = 1:numAng
        proj    = sinogram(:, angIdx).';
        projFFT = fft(proj, N_pad);
        projFFT = projFFT .* filterCurve;
        projF   = real(ifft(projFFT, N_pad));
        filteredSino(:, angIdx) = projF(1:numDet).';
    end
end

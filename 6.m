clear; clc; close all;

angleSets = {
    0:1:180,  ... 
    0:1:120,  ...
    0:1:90,   ...
    0:1:60,   ...
    0:1:30
};

angleLabels = { ...
    'Full [0,180°]', ...
    'Angular range [0,120°]', ...
    'Angular range [0,90°]', ...
    'Angular range [0,60°]', ...
    'Angular range [0,30°]'
};

N = 256;
phantomImg = phantom(N);

figure('Name','Original Phantom');
imshow(phantomImg, []);
title('2D Phantom (Shepp-Logan)');
drawnow;

for i = 1:length(angleSets)
    angles_i = angleSets{i};

    sinogram_i = radon(phantomImg, angles_i);

    recon_i = FFBP(sinogram_i, angles_i);

    figure('Name',sprintf('Reconstruction: %s', angleLabels{i}));
    imshow(recon_i, []);
    title(sprintf('Reconstruction: %s', angleLabels{i}));
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Helper Function (if not already defined in your workspace)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function recon = FFBP(sinogram, theta)

    [numDet, numAng] = size(sinogram);
    N_pad = 2^nextpow2(2 * numDet);

    freqHalf = (0:(N_pad/2 - 1)) / N_pad;
    rampFilter = [freqHalf, fliplr(freqHalf)];

    filteredSino = zeros(size(sinogram));
    for k = 1:numAng
        proj    = sinogram(:,k).';
        projFFT = fft(proj, N_pad);
        filtFFT = projFFT .* rampFilter;
        projF   = real(ifft(filtFFT, N_pad));
        filteredSino(:,k) = projF(1:numDet).';
    end

    recon = iradon(filteredSino, theta, 'linear','none', 1.0, numDet);
end

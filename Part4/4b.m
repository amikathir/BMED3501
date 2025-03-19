clear; clc; close all;

N = 256; 
uniformValue = 1.0; % The uniform background intensity
phantomImg = ones(N) * uniformValue;

theta = 0:1:179;
[cleanSino, xp] = radon(phantomImg, theta);

rng(0);
noiseLevel = 0.02; % or any amplitude you like
noisySino = cleanSino + noiseLevel * max(cleanSino(:)) * randn(size(cleanSino));

[numDet, numAng] = size(noisySino);
N_pad = 2^nextpow2(2*numDet);


freqHalf = (0:(N_pad/2 - 1)) / N_pad;
rampHalf = freqHalf; % basic ramp

ramLak_full = [rampHalf, fliplr(rampHalf)];

cutoff = freqHalf(end);
sheppLogan_full = ramLak_full .* sincfreq( linspace(-0.5, 0.5, N_pad), cutoff );

hann_full = ramLak_full .* hannWindow( linspace(-0.5, 0.5, N_pad), cutoff );

filters = {ramLak_full, sheppLogan_full, hann_full};
filterNames = {'Ram-Lak','Shepp-Logan','Hann'};

reconResults = cell(size(filters));
for fIdx = 1:length(filters)
    currentFilter = filters{fIdx};
    filteredSino = zeros(size(noisySino));

    for angIdx = 1:numAng
        proj = noisySino(:,angIdx).';
        projFFT = fft(proj, N_pad);
        
        % Multiply by current filter
        filteredFFT = projFFT .* currentFilter;
        
        % Inverse FFT and truncate
        filteredProj = real(ifft(filteredFFT, N_pad));
        filteredSino(:, angIdx) = filteredProj(1:numDet).';
    end

    % Backproject with 'none' so we only use our custom filter
    recon = iradon(filteredSino, theta, 'linear','none', 1.0, N);
    reconResults{fIdx} = recon;
end

figure('Name','Noisy Reconstructions with Various Filters');
for fIdx = 1:length(filters)
    subplot(1, length(filters), fIdx);
    imshow(reconResults{fIdx}, []);
    title(sprintf('%s Filter', filterNames{fIdx}));
end
sgtitle('Comparison of Noisy Phantom Reconstructions');

roiSize = 100; % size of the central ROI
center = N/2;  
roiIndices = (center - roiSize/2 + 1) : (center + roiSize/2);

figure;
for fIdx = 1:length(filters)
    reconImg = reconResults{fIdx};
    
    roi = reconImg(roiIndices, roiIndices);
    roiMean = mean(roi(:));
    roiNoiseOnly = roi - roiMean;
    
    F2 = fft2(roiNoiseOnly);
    NPS2D = abs(F2).^2 / (roiSize^2);  % Normalized power
    
    nps1D = radialAverage(NPS2D);
    
    freqStep = 1/roiSize;
    maxFreqIndex = floor(length(nps1D)/2); 
    freqAxis = (0:(maxFreqIndex-1)) * freqStep;  % We only plot 0..Nyquist
    
    subplot(1, length(filters), fIdx);
    plot(freqAxis, nps1D(1:maxFreqIndex), 'LineWidth',1.3);
    grid on; xlabel('Spatial Frequency (cycles/pixel)'); ylabel('NPS');
    title(filterNames{fIdx});
    axis tight;
end
sgtitle('Noise Power Spectrum Comparison');

freqSym = linspace(-0.5, 0.5, N_pad); 
figure;
hold on;
for fIdx = 1:length(filters)
    plot(freqSym, abs(filters{fIdx}), 'LineWidth',1.5, 'DisplayName', filterNames{fIdx});
end
legend('Location','best');
xlabel('Normalized Frequency');
ylabel('Magnitude');
title('MTF Curves (Filter Frequency Responses)');
grid on;
hold off;

%% ---- HELPER FUNCTIONS ----
function y = sincfreq(freq, fcut)
% Shepp-Logan window factor = sinc(pi * freq / cutoff).
% Only apply within |freq| <= cutoff; zero otherwise.
    y = zeros(size(freq));
    mask = (abs(freq) <= abs(fcut));
    % scaled freq
    x = freq(mask) ./ fcut;
    % safe Sinc definition: sin(pi*x)/(pi*x)
    y(mask) = sin(pi*x) ./ (pi*x);
    % handle x=0 => sinc(0)=1
    y(abs(x) < 1e-12) = 1;
end

function y = hannWindow(freq, fcut)
    y = zeros(size(freq));
    mask = (abs(freq) <= abs(fcut));
    x = freq(mask) ./ fcut;
    y(mask) = 0.5*(1 + cos(pi*x));
end

function profile1D = radialAverage(image2D)
    [nr, nc] = size(image2D);
    cy = floor(nr/2)+1;
    cx = floor(nc/2)+1;
    maxRadius = floor(min(nr, nc)/2);
    
    sumR = zeros(maxRadius,1);
    countR = zeros(maxRadius,1);
    
    for r = 1:nr
        for c = 1:nc
            dy = r - cy;
            dx = c - cx;
            radius = sqrt(dx^2 + dy^2);
            bin = floor(radius) + 1;
            if bin <= maxRadius
                sumR(bin) = sumR(bin) + image2D(r,c);
                countR(bin) = countR(bin) + 1;
            end
        end
    end
    
    profile1D = sumR ./ max(countR, 1);
end

clear; clc; close all;

N = 256;                           
phantomImg = phantom(N);           
theta = 0:1:179;                   
[cleanSino, xp] = radon(phantomImg, theta);

rng(0);                           
noiseLevel = 0.02;
noisySino = cleanSino + noiseLevel * max(cleanSino(:)) * randn(size(cleanSino));

[numDet, numAng] = size(noisySino);
N_pad = 2^nextpow2(2*numDet);      

freqHalf = (0:(N_pad/2 - 1)) / N_pad;

ramLak_half = freqHalf;                          
ramLak_full = [ramLak_half, fliplr(ramLak_half)];

cutoff = freqHalf(end);
sheppLogan_full = ramLak_full .* sincfreq(linspace(-0.5,0.5,N_pad), cutoff);

hann_full = ramLak_full .* hannWindow(linspace(-0.5,0.5,N_pad), cutoff);

filters      = {ramLak_full, sheppLogan_full, hann_full};
filterNames  = {'Ram-Lak','Shepp-Logan','Hann'};

freqSym = linspace(-0.5, 0.5, N_pad);
figure('Name','Filter Frequency Responses');
hold on;
for i = 1:numel(filters)
    plot(freqSym, abs(filters{i}), 'LineWidth', 1.5, 'DisplayName', filterNames{i});
end
xlabel('Normalized Frequency');
ylabel('Magnitude');
title('Comparison of Filter Frequency Responses');
legend('Location','best');
grid on;
hold off;
drawnow;

reconResults = cell(size(filters));
for fIdx = 1:numel(filters)
    currentFilter = filters{fIdx};
    
    filteredSino = zeros(size(noisySino));
    
    for angIdx = 1:numAng
        proj    = noisySino(:,angIdx).';
        projFFT = fft(proj, N_pad);
        filtFFT = projFFT .* currentFilter;
        projF   = real(ifft(filtFFT, N_pad));
        filteredSino(:, angIdx) = projF(1:numDet).';
    end
    
    reconResults{fIdx} = iradon(filteredSino, theta, 'linear','none',1.0, N);
end

figure('Name','Reconstruction with Different Filters');
for i = 1:numel(filters)
    subplot(1, numel(filters), i);
    imshow(reconResults{i}, []);
    title(filterNames{i});
end
sgtitle('Noisy Phantom Reconstructions (Ramp vs. Windowed Filters)');
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = sincfreq(freq, fcut)
    y = zeros(size(freq));
    mask = (abs(freq) <= abs(fcut));
    x = freq(mask) ./ fcut;
    y(mask) = sin(pi*x) ./ (pi*x); 
    y(abs(x) < 1e-12) = 1;
end

function y = hannWindow(freq, fcut)
    y = zeros(size(freq));
    mask = (abs(freq) <= abs(fcut));
    x = freq(mask) ./ fcut;
    y(mask) = 0.5*(1 + cos(pi*x));
end

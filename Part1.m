clear; clc; close all;

N = 256;                 
N_pad = 2^nextpow2(2*N);
freq = (-N_pad/2 : N_pad/2 - 1)/N_pad;  

idealRamp = abs(freq);

halfRamp = (0 : (N_pad/2 - 1)) / N_pad;
discreteRamp = [halfRamp, fliplr(halfRamp)];

figure('Name','Ramp Filter Frequency Response');
plot(freq, idealRamp, 'LineWidth',1.5); hold on;
plot(freq, discreteRamp, '--','LineWidth',1.5);
xlabel('Normalized Frequency');
ylabel('Magnitude');
legend('Ideal |S|','Discrete Ramp Filter','Location','best');
title('Ramp Filter Frequency Response');
grid on;
drawnow;

s = linspace(-10,10,N);
projectionRaw = exp(-0.5*(s/2).^2);  

proj_fft = fftshift(fft(projectionRaw, N_pad));

rampFilter = discreteRamp;

filtered_fft = proj_fft .* rampFilter;

filtered_projection_full = ifft(ifftshift(filtered_fft), N_pad);
filtered_projection = real(filtered_projection_full(1:N));

figure('Name','Raw Projection');
plot(s, projectionRaw, 'LineWidth',1.5);
title('Raw Projection');
xlabel('s');
ylabel('g(s)');
drawnow;

figure('Name','Fourier Transform (Magnitude)');
freqAxis = fftshift(freq); % for consistent plotting
plot(freqAxis, abs(proj_fft), 'LineWidth',1.5);
title('Fourier Transform (Magnitude)');
xlabel('Normalized Frequency'); 
ylabel('|FFT|');
drawnow;

figure('Name','After Ramp Filtering (Magnitude)');
plot(freqAxis, abs(filtered_fft), 'LineWidth',1.5);
title('After Ramp Filtering (Magnitude)');
xlabel('Normalized Frequency'); 
ylabel('|Filtered FFT|');
drawnow;

figure('Name','Filtered Projection (Spatial Domain)');
plot(s, filtered_projection, 'LineWidth',1.5);
title('Filtered Projection (Spatial Domain)');
xlabel('s');
ylabel('g_{filtered}(s)');
drawnow;

phantomSize = 128;
phantomImg = phantom(phantomSize);

angles = 0:10:170;
[sinogram, xp] = radon(phantomImg, angles);

[numDet, numAng] = size(sinogram);
N_pad_sino = 2^nextpow2(2*numDet);
fftMagArray = zeros(N_pad_sino, numAng);

freqAxis2 = linspace(-0.5, 0.5, N_pad_sino);

for aIdx = 1:numAng
    proj = sinogram(:, aIdx);
    projFFT = fftshift(fft(proj, N_pad_sino));
    fftMagArray(:, aIdx) = abs(projFFT);
end

figure('Name','FFT Magnitude vs. Angle & Frequency');
imagesc(angles, freqAxis2, fftMagArray);
colormap jet; colorbar;
xlabel('Projection Angle (degrees)');
ylabel('Normalized Frequency');
title('Magnitude of 1D FFT of Each Projection');
axis xy;
drawnow;

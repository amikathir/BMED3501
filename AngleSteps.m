clear; clc; close all;

N = 256; % Phantom/image size
phantomImg = phantom(N);

F2_phantom = fftshift(fft2(phantomImg));

freqAxis = linspace(-0.5, 0.5, N);
[Xf, Yf] = meshgrid(freqAxis, freqAxis);

figure('Name','Phantom and 2D Fourier Magnitude');
subplot(1,2,1);
imshow(phantomImg, []);
title('Original Phantom');
subplot(1,2,2);
imshow(log(abs(F2_phantom)+1), []);
title('2D FFT of Phantom (log scale)');
drawnow;


theta0 = 45;
[sinogramSingle, xp] = radon(phantomImg, theta0);

projection = sinogramSingle(:,1);
nDet = length(projection);

projFFT = fftshift(fft(projection));
freqProj = linspace(-0.5, 0.5, nDet);

figure('Name','Central Slice Theorem');
imshow(log(abs(F2_phantom)+1), []);
title(sprintf('Central Slice at \\theta = %d^{\\circ}', theta0));
hold on;

r = linspace(-0.5, 0.5, nDet);  % radial freq axis
xf_line = r * cosd(theta0);
yf_line = r * sind(theta0);

xf_px = (xf_line + 0.5) * N; 
yf_px = (yf_line + 0.5) * N;
plot(xf_px, yf_px, 'r', 'LineWidth',1.5);

figure('Name','1D FFT of Projection');
plot(freqProj, abs(projFFT), 'LineWidth',1.5);
title('Magnitude of 1D FFT of Projection');
xlabel('Normalized Frequency'); ylabel('|FFT|');
grid on;
drawnow;


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
    sino = radon(phantomImg, angles);
    [numDet, numAng] = size(sino);
    
    subplot(1, numPlots, i);
    imagesc(log(abs(F2_phantom)+1)); axis image;
    hold on; title(sprintf('\\Delta\\theta = %d^{\\circ}', angles(2)-angles(1)));
    
    for k = 1:numAng
        theta_k = angles(k);
        freq_s  = linspace(-0.5, 0.5, numDet);
        xf = freq_s .* cosd(theta_k);
        yf = freq_s .* sind(theta_k);
        
        % Convert freq coords to pixel coords
        xf_px = (xf + 0.5) * N;
        yf_px = (yf + 0.5) * N;
        plot(xf_px, yf_px, 'w'); 
    end
    
    axis off;
end
sgtitle('Comparison of Fourier Coverage for Different Theta Increments');
drawnow;

anglesSet = {
    0:1:179,  ... 
    0:2:179,  ...  
    0:5:179,  ...  
    0:10:179  ...  
};

N = 256;
phantomImg = phantom(N);

figure('Name','Aliasing Artifacts with Zoomed View');
numSets = length(anglesSet);

roiSize = 50;                   
centerCoord = N/2;              
roiStart = round(centerCoord - roiSize/2) + 1;
roiEnd   = roiStart + roiSize - 1;

for i = 1:numSets
    angles = anglesSet{i};
    samplingStep = angles(2) - angles(1);
    
    sino = radon(phantomImg, angles);
    
    recon = iradon(sino, angles, 'linear','Ram-Lak', 1.0, N);

    subplot(2, numSets, i);
    imshow(recon, []);
    title(sprintf('Full View (%d\\circ steps)', samplingStep));

    subplot(2, numSets, i + numSets);
    zoomedRegion = recon(roiStart:roiEnd, roiStart:roiEnd);
    imshow(zoomedRegion, []);
    title(sprintf('Zoomed (%d\\circ)', samplingStep));
end

sgtitle('Aliasing Artifacts: Angular Undersampling (Full vs. Zoomed Views)');
drawnow;


anglesSet = {
    0:1:179,  ...  
    0:2:179,  ...  
    0:5:179,  ...  
    0:10:179  ... 
};

N = 256;
phantomImg = phantom(N);

roiSize = 50;                   
centerCoord = N/2;             
roiStart = round(centerCoord - roiSize/2) + 1;
roiEnd   = roiStart + roiSize - 1;

for i = 1:length(anglesSet)
    angles = anglesSet{i};
    samplingStep = angles(2) - angles(1);
    
    sino = radon(phantomImg, angles);
    
    recon = iradon(sino, angles, 'linear','Ram-Lak', 1.0, N);
    
    zoomedRegion = recon(roiStart:roiEnd, roiStart:roiEnd);
    
    figure('Name', sprintf('Full Reconstruction: %d-degree steps', samplingStep));
    imshow(recon, []);
    title(sprintf('Full View (%d\\circ steps)', samplingStep));
    drawnow;
    
    figure('Name', sprintf('Zoomed Reconstruction: %d-degree steps', samplingStep));
    imshow(zoomedRegion, []);
    title(sprintf('Zoomed View (%d\\circ steps)', samplingStep));
    drawnow;
end


% Frequency Analysis

clear;
close all;

%% Vertebra CT

% Load image
volume = double(loadCT());
sz = size(volume);
sz2 = round(sz / 2);

% display
figure; mprov(volume);
I = squeeze(volume(sz2(1), :, :));
figure; imagesc(I); colormap gray; colorbar; axis image;

% FFT
F = abs(fftshift(fft2(I)));
figure; imagesc(log(1+F)); colormap gray; colorbar; axis image;
figure; hold on;
plot(F(sz2(2), :));

% DWT
[LoD, HiD] = wfilters('haar', 'd');
[cA, cH, cV, cD] = dwt2(I, LoD, HiD, 'mode', 'symh');
figure;
subplot(2, 2, 1); imagesc(cA); colormap gray; colorbar; axis image;
subplot(2, 2, 2); imagesc(cH); colormap gray; colorbar; axis image;
subplot(2, 2, 3); imagesc(cV); colormap gray; colorbar; axis image;
subplot(2, 2, 4); imagesc(cD); colormap gray; colorbar; axis image;

%% MRI




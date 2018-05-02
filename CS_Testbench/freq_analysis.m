% Frequency Analysis

clear;
close all;

num_images = 37;
C = 0.1:0.1:0.9;
SSIMs = zeros(num_images, length(C));

tic;
for k = 1:num_images

    % Load image
    if k == 1 % vertebra CT
        volume = double(loadCT());
    elseif k == 2 % head CT
        v = load('data\Head CT\headCT.mat');
        volume = v.u;
    else
        v = load(['data\MRI\', num2str(k - 2), '_scan.mat']);
        volume = v.img;
    end
    
    sz = size(volume);
    sz2 = sz / 2;
%     figure; mprov(volume);

    % DFT of image
    F = fftshift(fftn(volume));
%     figure; mprov(log(1 + abs(F)));

    % Nyquist cutoffs to test
    N = round(C' * sz);
    [x, y, z] = meshgrid(1:sz(2), 1:sz(1), 1:sz(3));
    x = x - sz2(2);
    y = y - sz2(1);
    z = z - sz2(3);

    for i = 1:length(N)
        cutoff = (x / N(i, 2)).^2 + (y / N(i, 1)).^2 + (z / N(i, 3)).^2 <= 1;
        F_cutoff = F .* cutoff;
        I = abs(ifftn(ifftshift(F_cutoff)));
        SSIMs(k, i) = ssim(I, volume);
    end
end
toc;

% results
C = [0, C, 1];
SSIMs = horzcat(zeros(num_images, 1), SSIMs, ones(num_images, 1));
save('results\freq_analysis.mat', 'SSIMs', 'C');

idx = sort(randperm(35, 3));
figure; hold on;
plot(C, SSIMs(1, :), 'LineWidth', 2);
plot(C, SSIMs(2, :), 'LineWidth', 2);
plot(C, SSIMs(idx(1), :), 'LineWidth', 2);
plot(C, SSIMs(idx(2), :), 'LineWidth', 2);
plot(C, SSIMs(idx(3), :), 'LineWidth', 2);
legend('Vertebra CT', 'Head CT', ['Head MRI ', num2str(idx(1))], ...
    ['Head MRI ', num2str(idx(2))], ['Head MRI ', num2str(idx(3))], 'location', 'best');
xlabel('Nyquist cutoff');
ylabel('SSIM');
title('Fourier Frequency sparsity of images');
set(gca, 'fontsize', 16);
saveas(gcf, 'results\freq_analysis.png');

figure;
subplot(2,3,1); imagesc(flipud(squeeze(volume(round(sz2(1)), :, :)))); colormap gray; colorbar; axis image; axis off; title('Sagittal');
subplot(2,3,2); imagesc(squeeze(volume(:, round(sz2(2)), :))); colormap gray; colorbar; axis image; axis off; title('Coronal');
subplot(2,3,3); imagesc(squeeze(volume(:, :, round(sz2(3))))); colormap gray; colorbar; axis image; axis off; title('Axial');
subplot(2,3,4); imagesc(log(1+abs(squeeze(F(round(sz2(1)), :, :)))));
    colormap gray; colorbar; axis image; xlabel('f_y'); ylabel('f_z'); set(gca, 'XTick', [], 'YTick', []);
subplot(2,3,5); imagesc(log(1+abs(squeeze(F(:, round(sz2(2)), :)))));
    colormap gray; colorbar; axis image; xlabel('f_x'); ylabel('f_z'); set(gca, 'XTick', [], 'YTick', []);
subplot(2,3,6); imagesc(log(1+abs(squeeze(F(:, :, round(sz2(3)))))));
    colormap gray; colorbar; axis image; xlabel('f_x'); ylabel('f_y'); set(gca, 'XTick', [], 'YTick', []);
saveas(gcf, 'results\freq_MRI.png');

v = load('data\Head CT\headCT.mat');
volume = v.u;
F = fftshift(fftn(volume));
figure;
subplot(2,3,1); imagesc(flipud(squeeze(volume(256, :, :)))); colormap gray; colorbar; axis image; axis off; title('Sagittal');
subplot(2,3,2); imagesc(squeeze(volume(:, 256, :))); colormap gray; colorbar; axis image; axis off; title('Coronal');
subplot(2,3,3); imagesc(squeeze(volume(:, :, 90))); colormap gray; colorbar; axis image; axis off; title('Axial');
subplot(2,3,4); imagesc(log(1+abs(squeeze(F(256, :, :)))));
    colormap gray; colorbar; axis image; xlabel('f_y'); ylabel('f_z'); set(gca, 'XTick', [], 'YTick', []);
subplot(2,3,5); imagesc(log(1+abs(squeeze(F(:, 256, :)))));
    colormap gray; colorbar; axis image; xlabel('f_x'); ylabel('f_z'); set(gca, 'XTick', [], 'YTick', []);
subplot(2,3,6); imagesc(log(1+abs(squeeze(F(:, :, 90)))));
    colormap gray; colorbar; axis image; xlabel('f_x'); ylabel('f_y'); set(gca, 'XTick', [], 'YTick', []);
saveas(gcf, 'results\freq_headCT.png');


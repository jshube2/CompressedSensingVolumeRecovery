% Frequency Analysis

clear;
close all;

%% Vertebra CT

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

C = [0, C, 1];
SSIMs = horzcat(zeros(num_images, 1), SSIMs, ones(num_images, 1));
save('results\freq_analysis.mat', 'SSIMs', 'C');

figure; hold on;
plot(C, SSIMs(1, :), 'LineWidth', 2);
plot(C, SSIMs(2, :), 'LineWidth', 2);
plot(C, SSIMs(3, :), 'LineWidth', 2);
plot(C, SSIMs(15, :), 'LineWidth', 2);
plot(C, SSIMs(30, :), 'LineWidth', 2);
legend('Vertebra CT', 'Head CT', 'Head MRI 1', 'Head MRI 13', 'Head MRI 28', 'location', 'southeast');
xlabel('Nyquist cutoff');
ylabel('SSIM');
title('Fourier Frequency sparsity of images');
set(gca, 'fontsize', 16);

%% MRI




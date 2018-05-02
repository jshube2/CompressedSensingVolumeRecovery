clear;
close all;
clc;
%-------------------------------------------------------------------------%
%-------------------------EXPERIMENT PARAMETERS---------------------------%
%-------------------------------------------------------------------------%
sparsity = [10 20 30 40 50 60 70 80 90];
%sparsity = 50;
modality = 'MRI';
recovery = 'Sparse';

%-------------------------------------------------------------------------%
%-------------------------RANDOM VOXEL DELETION---------------------------%
%-------------------------------------------------------------------------%

if strcmp(modality, 'CT')
    volume = loadCT();
end

if strcmp(modality, 'MRI')
    load('data/MRI/1_scan.mat');
    %volume = double(img);
    r=1;
    for i=1:size(img,3)
       if sum(img(:,:,i))==0
           fprintf('%g\n',i);
           continue;
       end
       volume(:,:,r) = img(:,:,i);
       r=r+1;
    end
end

figure; mprov(volume); % viewCrossSection(volume);
title('Original Volume');
%-------------------------------------------------------------------------%
%------------------------------(ALGORITHM A)------------------------------%
%-----------------------------------ALM-----------------------------------%
for i = 1:length(sparsity)
    switch recovery
        case 'Corrupt'
            sparse_volume = makeCorrupt(volume, sparsity(i)); 
        case 'Sparse'
            sparse_volume = makeSparse(volume, sparsity(i)); 
        case 'Patch'
            sparse_volume = makeSparsePatches(volume, sparsity(i), 4); 
    end
    
    %viewCrossSection(sparse_volume);
    %title('Sparse Volume');
    %makeVolumeMovie(sparse_volume); % USE THIS TO SAVE 3D VOLUME TO VIDEO
    
    %-DO-RECOVERY-HERE-%
    [recovered_volume, recovered_error] = recoveryAlgorithm(double(sparse_volume),'ALM',0.1);
    %recovered_volume = sparse_volume;
    %viewCrossSection(recovered_volume);
    %title('Recovered Volume');
    
    %------------------%
    error(i) = ssim(double(volume), recovered_volume);
end

figure;
plot(sparsity,error,'o-');
title(sprintf('ALM MC Recovery - %s',recovery));
xlabel('Percentage of Missing Data');
ylabel('SSIM');
xlim([0 100]);

fname = sprintf('Results/%s_Error_%s_MC',modality,recovery);
saveas(gcf, [fname '.png']);
save([fname '.mat'],'error');

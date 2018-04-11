clear all;
close all;
clc;
%-------------------------------------------------------------------------%
%-------------------------EXPERIMENT PARAMETERS---------------------------%
%-------------------------------------------------------------------------%
sparsity = [10 20 30 40 50 60 70 80 90];


%-------------------------------------------------------------------------%
%-------------------------RANDOM VOXEL DELETION---------------------------%
%-------------------------------------------------------------------------%
volume = loadCT();
viewCrossSection(volume);
%-------------------------------------------------------------------------%
%------------------------------(ALGORITHM A)------------------------------%
%-------------------------------------------------------------------------%
for i = 1:length(sparsity)
    sparse_volume = makeSparse(volume, sparsity(i));
    viewCrossSection(sparse_volume);
    %makeVolumeMovie(sparse_volume); % USE THIS TO SAVE 3D VOLUME TO VIDEO
    %-DO-RECOVERY-HERE-%
    recovered_volume = sparse_volume;
    %------------------%
    error(i) = immse(volume, recovered_volume);
end

plot(sparsity,error,'o-');
title('Algorithm A');
xlabel('Percentage of Missing Data');
ylabel('Mean Squared Error');
xlim([0 100]);

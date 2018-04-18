clear all;
close all;
clc;
%-------------------------------------------------------------------------%
%-------------------------EXPERIMENT PARAMETERS---------------------------%
%-------------------------------------------------------------------------%
%sparsity = [10 20 30 40 50 60 70 80 90];
sparsity = [50];


%-------------------------------------------------------------------------%
%-------------------------RANDOM VOXEL DELETION---------------------------%
%-------------------------------------------------------------------------%
volume = loadCT();
viewCrossSection(volume);
title('Original Volume');
%-------------------------------------------------------------------------%
%------------------------------(ALGORITHM A)------------------------------%
%-----------------------------------ALM-----------------------------------%
for i = 1:length(sparsity)
    %sparse_volume = makeCorrupt(volume, sparsity(i));
    sparse_volume = makeSparse(volume, sparsity(i));
    %sparse_volume = makeSparsePatches(volume, sparsity(i), 4);
    viewCrossSection(sparse_volume);
    title('Sparse Volume');
    %makeVolumeMovie(sparse_volume); % USE THIS TO SAVE 3D VOLUME TO VIDEO
    %-DO-RECOVERY-HERE-%
    [recovered_volume, recovered_error] = recoveryAlgorithm(double(sparse_volume),'ALM',0.1);
    %recovered_volume = sparse_volume;
    viewCrossSection(recovered_volume);
    title('Recovered Volume');
    %------------------%
    error(i) = immse(double(volume), recovered_volume);
end

figure;
plot(sparsity,error,'o-');
title('Algorithm A');
xlabel('Percentage of Missing Data');
ylabel('Mean Squared Error');
xlim([0 100]);

%-------------------------------------------------------------------------%
%------------------------------(ALGORITHM B)------------------------------%
%----------------------------------SVT------------------------------------%
for i = 1:length(sparsity)
    sparse_volume = makeSparse(volume, sparsity(i));
    viewCrossSection(sparse_volume);
    title('Sparse Volume');
    %makeVolumeMovie(sparse_volume); % USE THIS TO SAVE 3D VOLUME TO VIDEO
    %-DO-RECOVERY-HERE-%
    [recovered_volume, recovered_error] = recoveryAlgorithm(double(sparse_volume),'SVT',0.1);
    %recovered_volume = sparse_volume;
    viewCrossSection(recovered_volume);
    title('Recovered Volume');
    %------------------%
    error(i) = immse(double(volume), recovered_volume);
end

figure;
plot(sparsity,error,'o-');
title('Algorithm B');
xlabel('Percentage of Missing Data');
ylabel('Mean Squared Error');
xlim([0 100]);

%-------------------------------------------------------------------------%
%------------------------------(ALGORITHM C)------------------------------%
%----------------------------------APG------------------------------------%
for i = 1:length(sparsity)
    sparse_volume = makeSparse(volume, sparsity(i));
    viewCrossSection(sparse_volume);
    title('Sparse Volume');
    %makeVolumeMovie(sparse_volume); % USE THIS TO SAVE 3D VOLUME TO VIDEO
    %-DO-RECOVERY-HERE-%
    [recovered_volume, recovered_error] = recoveryAlgorithm(double(sparse_volume),'APG',0.1);
    %recovered_volume = sparse_volume;
    viewCrossSection(recovered_volume);
    title('Recovered Volume');
    %------------------%
    error(i) = immse(double(volume), recovered_volume);
end

figure;
plot(sparsity,error,'o-');
title('Algorithm B');
xlabel('Percentage of Missing Data');
ylabel('Mean Squared Error');
xlim([0 100]);

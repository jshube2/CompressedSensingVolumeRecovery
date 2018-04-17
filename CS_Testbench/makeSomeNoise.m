function noisy_volume = makeSomeNoise(volume, corruption)
    % INPUT: 
    %   volume = 3D image volume
    %   corruption = percentage of variance of voxel values [0 100]
    % OUTPUT:
    %   noisy_volume = 3D image volume with voxels randomly noised
    corruption = (corruption / 100) * mean(volume(:));
    noisy_volume = double(volume) + corruption*randn(size(volume));
end
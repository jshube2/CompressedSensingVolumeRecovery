function [A_hat, E_hat] = recoveryAlgorithm(volume, option, lambda)
    % INPUT:
    %   volume = 3D volume
    %   option = algorithm choice (i.e., 'ALM'...)
    %   lambda = sparsity penality parameter
    % OUTPUT:
    %   A_hat = recovered 3D volume
    %   E_hat = recovered 3D error matrix
    
    tol = 1e-6;
    maxIter = 1000;
    
    [num_x, num_y, num_z] = size(volume);

    % Preallocation
    A_hat = zeros(size(volume));
    E_hat = zeros(size(volume));

    switch option
        case 'ALM'
            addpath([cd '/inexact_alm_rpca']);
            % Apply RPCA via inexact ALM on each matrix through depth, z
            for z = 1:num_z
               [A_hat(:,:,z), E_hat(:,:,z)] = inexact_alm_rpca(volume(:,:,z), lambda, tol, maxIter);
            end
            
        case 'ALM-MC'
            addpath([cd '/inexact_alm_mc']);
            % Apply RPCA via inexact ALM on each matrix through depth, z
            for z = 1:num_z
               [A_hat(:,:,z), E_hat(:,:,z)] = inexact_alm_mc(volume(:,:,z), tol, maxIter);
            end
    end
end
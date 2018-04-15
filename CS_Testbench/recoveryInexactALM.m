function [A_hat, E_hat] = recoveryInexactALM(volume, lambda)
    % INPUT:
    %   volume = 3D volume
    %   lambda = 
    % OUTPUT:
    %   A_hat = recovered 3D volume
    %   E_hat = recovered 3D error matrix
    
    addpath([cd '/inexact_alm_rpca']);
    
    tol = 1e-6;
    maxIter = 1000;
    
    [num_x, num_y, num_z] = size(volume);

    % Preallocation
    A_hat = zeros(size(volume));
    E_hat = zeros(size(volume));

    % Apply RPCA via inexact ALM on each matrix through depth, z
    for z = 1:num_z
       [A_hat(:,:,z), E_hat(:,:,z)] = inexact_alm_rpca(volume(:,:,z), lambda, tol, maxIter);
    end

end
function sparse_volume = makeSparse(volume, sparsity)
    % INPUT: 
    %   volume = 3D image volume
    %   sparsity = Percentage of voxels to be removed [0 100]
    % OUTPUT:
    %   sparse_volume = 3D image volume with voxels randomly deleted
    
    sparse_volume = volume;
    
    [y_size, x_size, z_size] = size(volume);
    num_voxels = y_size*x_size*z_size;
    
    p = sparsity / 100;
    
    num2rmv = p*num_voxels;

    for ii = 1:y_size
        for jj = 1:x_size
            for kk = 1:z_size
                x = rand(1);
                if x < p
                    sparse_volume(ii,jj,kk) = 0;
                    num2rmv = num2rmv - 1;
                    sprintf("Made element into zero");
                end
                num_voxels = num_voxels - 1;
                p = num2rmv / num_voxels;
            end
        end
    end
end
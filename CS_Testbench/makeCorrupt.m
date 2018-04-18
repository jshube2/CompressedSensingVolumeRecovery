function crrpt_volume = makeCorrupt(volume, corruption)
    % INPUT: 
    %   volume = 3D image volume
    %   corruption = percentage of voxels to corrupt [0 100]
    % OUTPUT:
    %   crrpt_volume = 3D image volume with voxels randomly corrupted

    crrpt_volume = volume;
    max_val = max(volume(:));
    
    [y_size, x_size, z_size] = size(volume);
    num_voxels = y_size*x_size*z_size;
    
    p = corruption / 100;
    
    num2rmv = p*num_voxels;

    for ii = 1:y_size
        for jj = 1:x_size
            for kk = 1:z_size
                x = rand(1);
                if x < p
                    crrpt_volume(ii,jj,kk) = max_val * rand(1);
                    num2rmv = num2rmv - 1;    
                end
                num_voxels = num_voxels - 1;
                p = num2rmv / num_voxels;
            end
        end
    end
end
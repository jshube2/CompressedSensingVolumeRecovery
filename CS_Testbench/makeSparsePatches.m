function sparse_volume = makeSparsePatches(volume, sparsity, patch_size)
    % INPUT: 
    %   volume = 3D image volume
    %   sparsity = Percentage of voxels to be removed [0 100]
    %   patch_size = edge length of n x n x n cube to be deleted
    % OUTPUT:
    %   sparse_volume = 3D image volume with voxel patches randomly deleted
    
    sparse_volume = volume;
    
    [y_size, x_size, z_size] = size(volume);
    num_voxels = y_size*x_size*z_size;
    
    p = sparsity / 100;
    cube_size = patch_size ^ 3;
    num2rmv = p*num_voxels/cube_size;

    for ii = 1 + int16(patch_size/2):y_size - int16(patch_size/2)
        for jj = 1 + int16(patch_size/2):x_size - int16(patch_size/2)
            for kk = 1 + int16(patch_size/2):z_size - int16(patch_size/2)
                x = rand(1);
                if x < p && neighborsValid(sparse_volume, ii, jj, kk, patch_size)
                    % Delete patch
                    for iii = ii-int16(patch_size/2):ii+int16(patch_size/2)
                        for jjj = jj-int16(patch_size/2):jj+int16(patch_size/2)
                            for kkk = kk-int16(patch_size/2):kk+int16(patch_size/2)
                                sparse_volume(iii,jjj,kkk) = 0;
                            end
                        end
                    end
                    %
                    num2rmv = num2rmv - 1;
                end
                num_voxels = num_voxels - 1;
                p = num2rmv / (num_voxels);
            end
        end
    end
end

function result = neighborsValid(sparse_volume, ii, jj, kk, patch_size)
    for iii = ii-int16(patch_size/2):ii+int16(patch_size/2)
        for jjj = jj-int16(patch_size/2):jj+int16(patch_size/2)
            for kkk = kk-int16(patch_size/2):kk+int16(patch_size/2)
                if sparse_volume(iii,jjj,kkk) == 0
                    result = false;
                    return
                end
            end
        end
    end
    result = true;
    return
end

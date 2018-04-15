function volume = loadCT()
    % INPUT: 
    %   N/A
    % OUTPUT:
    %   volume = 3D CT volume
    
    rootdir = pwd;
    % Read CT Data
    cd([rootdir '/data/CT']);
    volumeList = dir('*IMA');
    for i = 1:length(volumeList)
        volume_raw(:,:,i) = dicomread(volumeList(i).name);
    end
    cd(rootdir);
    volume = volume_raw(127:256, :, :);
end
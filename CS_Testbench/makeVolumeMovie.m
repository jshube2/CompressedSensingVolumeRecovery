function makeVolumeMovie(volume, xsf, ysf, zsf, filename)
    % INPUT: 
    %   volume = 3D image volume
    %   xsf = scale factor in x direction
    %   ysf = scale factor in y direction
    %   zsf = scale factor in z direction
    %   filename = name of video to save
    % OUTPUT:
    %   Video file saved to location 'filename'
    
    if nargin < 2
        xsf = 1/0.00033333; %[px / m]
        ysf = 1/0.00033333;
        zsf = 1/0.00033333;
        filename = 'volume.mp4';
    end

    [y_size, x_size, z_size] = size(volume);

    x = 100*[1:x_size]/xsf;
    y = 100*[1:y_size]/ysf;
    z = 100*[1:z_size]/zsf;
    [yg, xg, zg] = meshgrid(x,y,z);

    %Make video
    f1 = figure;
    set(gcf,'color','w');
    set(f1,'Position',[100,100,800,800]);
    vidfile = VideoWriter([pwd '\vid\' filename]);
    vidfile.FrameRate = 20;
    open(vidfile);

    for i = 1:1:y_size
        slice(yg,xg,zg,double(volume), [], [y(i)], [0])
        shading interp
        colormap gray
        view([-144 40])
        axis image
        xlabel('Lateral [cm]')
        ylabel('Axial [cm]')
        zlabel('Elevation [cm]')
        title('Volume')
        drawnow

        writeVideo(vidfile,getframe(f1));
        figure(f1);clf
    end

    close(vidfile)
end
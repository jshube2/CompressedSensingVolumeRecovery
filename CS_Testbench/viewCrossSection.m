function viewCrossSection(volume, xsf, ysf, zsf)
    % INPUT: 
    %   volume = 3D image volume
    %   xsf = scale factor in x direction
    %   ysf = scale factor in y direction
    %   zsf = scale factor in z direction
    % OUTPUT:
    %   Simple view of the 3D volume with 3 plane slice view
    
    if nargin < 2
        xsf = 1/0.00033333; %[px / m]
        ysf = 1/0.00033333;
        zsf = 1/0.00033333;
    end

    [y_size, x_size, z_size] = size(volume);

    x = 100*[1:x_size]/xsf;
    y = 100*[1:y_size]/ysf;
    z = 100*[1:z_size]/zsf;
    [yg, xg, zg] = meshgrid(x,y,z);
    figure
    slice(yg,xg,zg,double(volume), [x(x_size/2)], [y(y_size/2)], [z(z_size/2)])
    shading interp
    colormap gray
    view([-144 40])
    axis image
    xlabel('Lateral [cm]')
    ylabel('Axial [cm]')
    zlabel('Elevation [cm]')
    title('Volume')
    drawnow
end
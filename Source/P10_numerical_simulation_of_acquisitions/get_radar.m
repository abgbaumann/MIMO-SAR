function [X_North,Y_East,Z_Height,...
          x_along_axis,y_cross_axis,...
          theta_c,theta_l,theta_r,...
          range_c,range_u,range_d] = get_radar(X0,Y0,Z0,Azi)

    % This functions returns the SLC grid of simulated radar data based on
    % the instruments coordinates (X0,Y0,Z0) and its Azimuth.
    %
    % -------------------------------------------------------------------------
    % by Andreas Baumann-Ouyang, ETH Zürich, GSEG (8th September 2022)
    %

    plt_fig = 0;
    
    numAngles = 128;
    numRange = 512;
    rngRes = 0.04;

    theta = asin(-2*((-numAngles/2:numAngles/2)/numAngles));
    dtheta = diff(theta);
    theta_l = theta(1:end-1);
    theta_r = theta(1:end-1) + dtheta;
    theta_c = (theta_l + theta_r) / 2;
    
    range_c = double(1:numRange) * rngRes - rngRes/2;
    range_d = range_c - rngRes/2;
    range_u = range_c + rngRes/2;
    
    [range_mat, theta_mat] = meshgrid(range_c,theta_c);
    
    x_along_axis = range_mat.*cos(theta_mat);
    y_cross_axis = range_mat.*sin(theta_mat);

    R = [cos(Azi*pi/180), sin(Azi*pi/180);...
         -sin(Azi*pi/180),cos(Azi*pi/180)];

    y_cross_axis = y_cross_axis(:);
    x_along_axis = x_along_axis(:);
        
    C_loc = [y_cross_axis,x_along_axis]';

    C_glo = R * C_loc + [Y0;X0]; 

    Y_East = C_glo(1,:);
    X_North = C_glo(2,:);
    Z_Height = ones(size(Y_East)) * Z0;
    
    if plt_fig
        figure;
        scatter(y_cross_axis(:),x_along_axis(:),'or',...
            'filled','DisplayName','Original');
        hold on
        axis equal
        scatter(Y_East,X_North,'ok',...
            'filled','DisplayName','Transformed');
        legend('location','best');
        box on
        grid on
    end

end


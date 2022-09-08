function [X,Y,Z] = gaussian_spatial_3d_object(y0,y1,dy,z0,z1,dz,xmax,N,X0,plt_fig)
% This functions simulates a deformation along a linear object (e.g.
% bridge).
%
% -------------------------------------------------------------------------
%
% y0:           [float]     - Lower bound of the grid on the Y-Axis
% y1:           [float]     - Upper bound of the grid on the Y-Axis
% dy:           [float]     - Step sizes of the grid on the Y-Axis
% z0:           [float]     - Lower Bound of the grid on the Z-Axis
% z1:           [float]     - Upper Bound of the grid on the Z-Axis
% dz:           [float]     - Step sizes of the grid on the Z-Axis
% xmax:         [float]     - Maximal anticipated displacement
% N:            [integer]   - Number of simulations
% X0:           [float]     - Average position of the grid on the X-Axis
% plt_fig       [0] or [1]  - Define if figures should be plotted
%
% -------------------------------------------------------------------------
%
% X,Y,Z:        [p x t]     - Time series of positions for each of the axes
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich, GSEG (8th September 2022)
%

Y = y0:dy:y1;
N_y = length(Y);
Z = z0:dz:z1;
N_z = length(Z);

% Spatial
Gauss_fact1 = 3;
ddy = (max(Y)-min(Y))/Gauss_fact1;
X = gaussmf(Y,[ddy median(Y)]); % get gaussian distribution
X = X - min(X); % Shift ends to zero
X = X / max(X) * xmax; % Scale to max displacement

% Temporal
Up_factor = 0.2;
Gauss_fact2 = 6;
N_third = floor(N/3);
N_center = N - (2*N_third);
T_third = [0:N_third-1];
dTup = (max(T_third)-min(T_third))/Gauss_fact2;
Tmed = median(T_third);
T_up = gaussmf(T_third,[dTup Tmed]); % Gaussian distribution (Up)
T_up = T_up - min(T_up);
T_up = T_up / max(T_up) * Up_factor; % Adjust for zero at ends and max
T_center = [0:N_center-1];
dTdown = (max(T_center)-min(T_center))/Gauss_fact2;
Tmed = median(T_center);
T_down = gaussmf(T_center,[dTdown Tmed]);  % Gaussian distribution (Down)
T_down = T_down - min(T_down);
T_down = T_down / max(T_down) * -1; % Adjust for zero at ends and max

T_scale = [T_up,T_down,T_up];
T_scale = T_scale / max(abs(T_scale));

% Combination
X = repmat(X,N,N_z) * -1;
X = X .* repmat(T_scale',1,N_y*N_z);

[Z,Y] = meshgrid(Z,Y);

X = round(X,5) + X0;
Z = repmat(Z(:),1,N)';
Y = repmat(Y(:),1,N)';

if plt_fig
    figure;
    hold on
    xlim([min(Y(:)),max(Y(:))]);
    ylim([min(X(:)),max(X(:))]);
    zlim([min(Z(:)),max(Z(:))]);
    view(-30,40);
    scatter3(Y(1,:),X(1,:),Z(1,:),'ok');
    for t_i = 1:N
        title(t_i)
        for pt_i = 1:N_y*N_z
            plt_obj(pt_i) = plot3([Y(1,pt_i),Y(t_i,pt_i)],...
                                    [X(1,pt_i),X(t_i,pt_i)],...
                                    [Z(1,pt_i),Z(t_i,pt_i)],...
                                    'r',...
                                    'LineWidth',3);
        end
        pause(0.1);
        delete(plt_obj);
    end
end

end


function [X,Y,Z] = sinusodial_time_3d_object(y0,y1,dy,z0,z1,dz,xmax,N,X0,lambda,N_waves,plt_fig)
% This functions simulates an oscillating.
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
% lambda:       [float]     - Damping factor
% N_waves:      [integer]   - Number of waves.
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

T = linspace(0,1,N);
Tw = T*2*pi; % "Time"

X = exp(-lambda.*Tw) .* sin(N_waves * Tw);
X = X / max(abs(X)) * xmax;

[Z,Y] = meshgrid(Z,Y);

X = repmat(X(:),1,N_y*N_z);
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


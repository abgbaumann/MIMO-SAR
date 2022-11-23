function [X,Y,Z] = planar_3d_object(y0,y1,dy,x,z0,z1,dz,xmax,N,X0,plt_fig,ax)
% This functions simulates a deformation defined by a parametric function.
%
% -------------------------------------------------------------------------
%
% y0:           [float]     - Lower bound of the grid on the Y-Axis
% y1:           [float]     - Upper bound of the grid on the Y-Axis
% dy:           [float]     - Step sizes of the grid on the Y-Axis
% x:            [sym]       - Parametric function (e.g. x = @(a, z) a*z;)
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
Z = z0:dz:z1;

A = linspace(0,xmax,N) ;
T = 0:N-1;

if ~exist('ax','var') || strcmp(ax,'Z')
    for z_i = 1:length(Z)
        for t_i = 1:length(A)
            X(t_i,:,z_i) = repmat(x(T(t_i),Z(z_i)),1,length(Y));
        end
    end
elseif strcmp(ax,'Y')
    for y_i = 1:length(Y)
        for t_i = 1:length(A)
            X(t_i,y_i,:) = repmat(x(T(t_i),Y(y_i)),1,length(Z));
        end
    end
else
    error('Ax unknown\n');
end

X_min = squeeze(min(X,[],1));
X_max = squeeze(max(X,[],1));
dX = abs(X_max - X_min);
scale_factor = xmax / max(dX(:));
X = X * scale_factor;
X = round(X,5) + X0;

if plt_fig
    figure;
    hold on
    xlim([min(Y),max(Y)]);
    ylim([min(X(:)),max(X(:))]);
    zlim([min(Z(:)),max(Z(:))]);
    view(-30,40);
    for t_i = 1:length(A)
        title(t_i)
        for z_i = 1:length(Z)
            plt_obj(z_i) = plot3(Y,X(t_i,:,z_i),repmat(Z(z_i),1,length(Y)),'r','LineWidth',5);
        end
        pause(0.1);
        delete(plt_obj);
        for z_i = 1:length(Z)
            plot3(Y,X(t_i,:,z_i),repmat(Z(z_i),1,length(Y)),'k');
        end
        pause(0.1);
    end
end

Y = repmat(Y,size(X,1),1,size(Z,2));
Z = permute(repmat(Z',1,size(X,1),size(X,2)),[2,3,1]);

X = reshape(X, [size(X,1), size(X,2)*size(X,3)]);
Y = reshape(Y, [size(X,1), size(X,2)*size(X,3)]);
Z = reshape(Z, [size(X,1), size(X,2)*size(X,3)]);

end


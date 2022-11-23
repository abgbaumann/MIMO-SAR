function [x_along,y_cross,z_height] = cart2cart_proj(X_north,Y_east,Z_height,...
                                                     X_inst,Y_inst,Z_inst,...
                                                     az_inst,el_inst,rot_inst,unit)
% This function takes a list of 3D coordinates and instrument properties 
% (coordinates/angles) in the same (Global) coordinate systems (CS) and 
% transforms it to local instrument coordinates.
%
% -------------------------------------------------------------------------
% INPUT:
% X_north, Y_east, Z_height:     Global Coordinates of the Point Cloud
% X_inst, Y_inst, Z_inst:        Global Cordinates of the Instrument
% az_inst, el_inst, rot_inst:    Rotation Angles of the Instrument
%
% OUTPUT:
% x_along, y_cross, z_height:    Point Cloud in the Local Instrument 
%                                Coordinate System.
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich, GSEG (12th March 2021)
%
    
% X_north = [0,2];
% Y_east = [0,0];
% Z_height = [0,0];
% 
% X_inst = 5;
% Y_inst = 0;
% Z_inst = 0;
% 
% az_inst = 180;
% el_inst = 45;
% rot_inst = 0;
% 
% unit = 'deg';

plt_fig = 0; % for ploting set plot_figure = 1

%% Translation:
% x/y/z is 0, if instrument CS is the same for the point cloud.
x_inst = 0;
y_inst = 0;
z_inst = 0;

X_trans = x_inst - X_inst;
Y_trans = y_inst - Y_inst;
Z_trans = z_inst - Z_inst;

T = eye(4);
T(1:3,4) = [X_trans; Y_trans; Z_trans]; % Translation Matrix

%% Rotation:
if strcmp(unit,'deg') || strcmp(unit,'degree')
    scale = pi / 180;
elseif strcmp(unit,'rad')
    scale = 1;
elseif strcmp(unit,'grad')
    scale = pi / 200;
end

alpha = rot_inst * scale;
beta = el_inst * scale;
gamma = az_inst * scale;

Rx = [1,          0,           0;...
      0, cos(alpha), -sin(alpha);...
      0, sin(alpha),  cos(alpha)];
  
Ry = [ cos(beta), 0,   sin(beta); ...
               0, 1,           0;...
      -sin(beta), 0,   cos(beta)];
  
Rz = [ cos(gamma), sin(gamma),   0;...
      -sin(gamma),  cos(gamma),   0;...
      0,                     0,  1]; % changed (clockwise)
  
R = eye(4);
R(1:3,1:3) = Rx * Ry * Rz; % Rotation Matrix

%% Coordinates:
X = [X_north; Y_east; Z_height; ones(size(Z_height))]; % Homogenous Coord.

%% Coordinate Transformation:
x_proj = R*T*X; % Transformation
x = x_proj(1:3,:) ./ x_proj(4,:); % Homogenous coordinates to 3D Coord.

%% Split
x_along = x(1,:);
y_cross = x(2,:);
z_height = x(3,:);

if plt_fig
    figure('units','normalized','outerposition',[0 0 1 1]);
    scatter3(y_cross(x_along>0),x_along(x_along>0),z_height(x_along>0),[],'b','filled');
    hold on
    scatter3(y_inst,x_inst,z_inst,50,'r','filled');
    xlabel('Cross [m]');
    ylabel('Along [m]');
    zlabel('Height [m]');
    view(0,0);
    pause(10);
    axis equal
end

end


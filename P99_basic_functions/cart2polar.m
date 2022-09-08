function [range, phi, theta] = cart2polar(x_along, y_cross, z_height, unit)

% This function takes a list of 3D coordinates and calculates the range,
% azimuth, and elevation angle with respect to the origin (0,0,0).
%
% -------------------------------------------------------------------------
% INPUT:
% x_along, y_cross, z_height:    Coordinates of the Point Cloud.
% unit:                          Units ('deg','rad', or 'grad') of the
%                                output angles.
%
% OUTPUT:
% range, phi, theta:             Range, Azimuth (phi), and Elevation 
%                                (theta) of the  Point Cloud with respect
%                                to the origin.
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich, GSEG (12th March 2021)
%

% unit='deg';
% x_along = rand(5,1)*20;
% y_cross = rand(5,1)*20;
% z_height = rand(5,1)*20;

range = sqrt(x_along.^2 + y_cross.^2 + z_height.^2);
phi = atan2(y_cross, x_along);
theta = sin(z_height ./ range);

if strcmp(unit,'deg') || strcmp(unit,'degree')
    scale = 180 / pi;
    digits = 8;
elseif strcmp(unit,'rad')
    scale = 1;
    digits = 10;
elseif strcmp(unit,'grad')
    scale = 200 / pi;
    digits = 8;
end
phi = round(phi * scale, digits);
theta = round(theta * scale, digits);
range = round(range,4);

end


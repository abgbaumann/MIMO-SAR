function [range_val, phi_val, ...
          range_c, phi_c, ...
          range_d, range_u, ...
          phi_l, phi_r,...
          y_cross_axis_c, x_along_axis_c] = cart2binrng2(y_cross_axis,...
                                                         x_along_axis,...
                                                         azimuths,...
                                                         ranges,...
                                                         unit_in,...
                                                         unit_out,....
                                                         numAngles)
% This function takes a vector of 2D coordiantes in y (cross/East-West) 
% and x (along/South-North) with given number of antennas and range 
% resolution and calculates the range and angle phi with
% the center position and left/right or up/down. The x/y coordinates
% correspond to the center positions of the bins.
%
% -------------------------------------------------------------------------
% INPUT:
% y_cross_axis, x_along_axis:    Vector of 2D coordinate of a point cloud
% azimuths:                      Vector of Azimuth Values.
% ranges:                        Vector of Range Values.
% unit:                          Indicating if the angles should be 
%                                returned in degrees (deg), radian (rad),
%                                or grad.
%
% OUTPUT:
% range_c:                       Range calculated from 0/0 coordinates.
%                                (Center)
% phi_c:                         Azimuth with phi = 0° is the direction of 
%                                Field of View and phi = 90° is orthogonal 
%                                to it. (Center)
% range_l,range_u:               Lower and upper boundary of the bin in 
%                                range.
% phi_l,phi_r:                   Left and right boundary of the bin in phi.
% y_cross_axis_cm, x_along_axis_c: Vector of 2D coordinate of a point cloud
%                                for the center coordinate of the bin.
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich, GSEG (14th June 2021)
%

plt_fig = 0;

phi_c = round(atan2(y_cross_axis,x_along_axis),6);

%% Angles / Azimuth
for i=1:2
    if i==1
        unit = unit_in;
    else
        unit = unit_out;
    end
    if strcmp(unit,'deg') || strcmp(unit,'degree')
        scale = 180 / pi;
        digits = 6;
    elseif strcmp(unit,'rad')
        scale = 1;
        digits = 10;
    elseif strcmp(unit,'grad')
        scale = 200 / pi;
        digits = 6;
    end
    if i==1
        scale_in = 1/scale;
    else
        scale_out = scale;
        digits_out = digits;
    end
end

if ~exist('numAngles','var')
    numAngles = length(azimuths);
end
phi_c_all = asin(-2*((-numAngles/2:numAngles/2)/numAngles));
dphi_all = diff(phi_c_all);
phi_l_all = phi_c_all(1:end-1);
phi_r_all = phi_c_all(1:end-1) + dphi_all;
phi_c_all = (phi_l_all + phi_r_all) / 2;

phi_c_all = round(phi_c_all,6);
dphi = zeros(length(phi_c),1);
for phi_all_i = 1:length(phi_c_all)
    [a,~] = ismember(phi_c,phi_c_all(phi_all_i));
    dphi(a) = dphi_all(phi_all_i);
end

phi_c = phi_c * scale_out;
dphi = dphi * scale_out;
phi_l = phi_c - dphi/2;
phi_r = phi_c + dphi/2;

phi_l = round(phi_l, digits);
phi_r = round(phi_r, digits);
phi_c = round(phi_c, digits);

phi_val = unique(phi_c);

%% Ranges
digits_dist = 3;
range_c = sqrt(x_along_axis.^2 + y_cross_axis.^2);
drange = median(ranges(2:end)-ranges(1:end-1));

range_d = range_c - drange/2;
range_u = range_c + drange/2;

range_c = round(range_c,digits_dist);
range_d = round(range_d,digits_dist);
range_u = round(range_u,digits_dist);

range_val = unique([range_u; range_d]);

y_cross_axis_c = range_c .* sin(phi_c*pi/180);
x_along_axis_c = range_c .* cos(phi_c*pi/180);

if plt_fig
    figure;
    ax = polaraxes;
    polarscatter(phi_c*pi/180,range_c,[],'^k'); hold on
    polarscatter(phi_l*pi/180,range_d,[],'ob'); hold on
    polarscatter(phi_r*pi/180,range_u,[],'xb'); hold on
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    ax.ThetaLim = [-90 90];
end

end



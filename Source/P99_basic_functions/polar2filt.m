function [idx] = polar2filt(range, phi, theta, drange, dphi, dtheta)
% This function takes a list of polar coordinates and computes which 
% points are within the (squared) field of view. 
%
% -------------------------------------------------------------------------
% INPUT:
% range, phi, theta:             Range, Azimuth (phi), and Elevation 
%                                (theta) of the  Point Cloud with respect
%                                to the origin.
% drange, dphi, dtheta:          Filtering values for Range (max.
%                                Distance), dphi (max. Azimuth), and dtheta 
%                                (max. Elevation).
%
% OUTPUT:
% idx:                           Indices of the points within the field of
%                                view.
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich, GSEG (12th March 2021)
%

idx = find(range >= min(drange) & range <= max(drange) & ...
           phi >= min(dphi) & phi <= max(dphi) & ...
           theta >= min(dtheta) & theta <= max(dtheta));

end


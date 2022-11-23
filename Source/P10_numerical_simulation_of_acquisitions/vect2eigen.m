function [amp, eigv] = vect2eigen(dXYZ)
% Takes a LOS-Displacement vector dXYZ and returns the eigen vector (unit
% vector) and signed absolute displacements.
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich, GSEG (8th September 2022)
%
Xpos = dXYZ(:,1) > 0;
Ypos = dXYZ(:,2) > 0;
Zpos = dXYZ(:,3) > 0;

dXYZFactor = rem(Xpos+Ypos+Zpos,2);
dXYZFactor(dXYZFactor==0) = -1;

amp = sqrt(sum(dXYZ.^2, 2)) .* dXYZFactor; % [L x 1] - Observed Displacement in LOS

eigv = dXYZ ./ amp; % Unit Vector of Displacement in 3D

eigv = round(eigv,5);

eigv = unique(eigv,'row');

eigv(logical(sum(isnan(eigv),2)/size(eigv,2)),:) = [];

if size(eigv,1)>1
    error('Attention, eigenvectors are not unique!')
elseif isempty(eigv)
    eigv = [nan nan nan];
end

end


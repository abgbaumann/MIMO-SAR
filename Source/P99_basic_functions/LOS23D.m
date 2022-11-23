function [x] = LOS23D(d_los, s_los)
% This functions takes as an input a coordinate vectors d_los with
% displacements in Line-of-Sight (LOS) and a uncertainty vector s_los 
% expressing the uncertainty of observations in LOS.
%
% -------------------------------------------------------------------------
%
% d_los:            [n x m] - Displacements as 3D vector
% s_los:            [n x 1] - Uncertainty (e.g. standard deviation) 
% x:                [1 x m] - Adjusted Coordinates
% 
% - with n being the number of observations (e.g. 3 different acquisitions)
% - with m being the dimension of the coordiante reference system (e.g. 3
%   for 3D)
%
% -------------------------------------------------------------------------
% 
% Coefficient Matrix:           A = d_los ./ L
% 
% Weigth Matrix:                P = diag( 1./ s_los.^2)
%
% Observation Vector:           L = sqrt( sum(d_los.^2, 2) )
%
% Observation Equation:         v = A * x - L
%
% Normal Equation:              N = A^T * P * A
%
% h-Vector:                     h = A^T * P * L
%
% Vector of Unknown Parameters: x = N^(-1) * h
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich, GSEG (07th March 2022)
%
  
plt_fig = 0;

L = sqrt(sum(d_los.^2, 2)); % [L x 1] - Observed Displacement in LOS

A = d_los ./ L; % Unit Vector of Displacement in 3D

% %% Adjusting of L for LSQ (after creating unit vector)
% XYZpos = sum(d_los./abs(d_los), 2);
% dXYZFactor = rem(XYZpos,2) * -1;
%         
% L = L .* dXYZFactor;

%% Test
if plt_fig
    figure; 
    hold on
    d1 = [zeros(size(d_los(:,1))),d_los(:,1)];
    d2 = [zeros(size(d_los(:,2))),d_los(:,2)];
    d3 = [zeros(size(d_los(:,3))),d_los(:,3)];
    for i=1:size(d_los(:,1),1)
        plot3(d1(i,:),d2(i,:),d3(i,:),'-k');
    end
    axis equal
end

P = diag(1./s_los.^2); % Uncertainty of the Observation (Matrix)

x = LSQ(L, A, P);


end


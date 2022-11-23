function [idx_azi,idx_rng,idx_azi_rng] = cart2bin(XR_N,YR_E,ZR_H,...
                                                  XO_N,YO_E,ZO_H,...
                                                  rng_bin_d,rng_bin_u,...
                                                  azi_bin_l,azi_bin_r,...
                                                  AZI_R,...
                                                  plt_val)
% This function takes a list of 3D coordinates of an instrument and a point
% cloud as well as polar coordiantes of the border points of all bins and
% returns the bin indice for each point in the point cloud
%
% -------------------------------------------------------------------------
% INPUT:
% XR_N,YR_E,ZR_H:                Global Coordinates of the Instrument
% XO_N,YO_E,ZO_H:                Global Cordinates of the Point Cloud
% rng_bin_d,rng_bin_u:           Polar Coordinate in Range for the lower
%                                (d) and upper (u) cornera of the bin.
% azi_bin_l,azi_bin_r:           Local Polar Coordinate in Azimuth for the 
%                                left (l) and right (r) corners of the bin 
%                                in [rad].
% AZI_R (optional):              Global Azimuth of the Instrument in [rad].
%
% OUTPUT:
% idx_azi,idx_rng,idx_azi_rng:   Indices in with respect to the bins.
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich, GSEG (23th March 2021)
%

    plt_fig = 0;
    scale_for_test = 0;
    
    if ~exist('AZI_R','var')
        AZI_R = 0;
    end
    
    dX = XO_N - XR_N; % Coord_Center - Coord_PointCloud
    dY = YO_E - YR_E;
    dZ = ZO_H - ZR_H;

    dRng = sqrt(dX.^2 + dY.^2 + dZ.^2);

    dAzi = atan2(dY,dX);

    dAzi = dAzi - AZI_R;% - AZI_R;
    
    N_theta = size(azi_bin_l,2);
    N_points = size(dAzi,2);

    if plt_fig
        if ~exist('plt_val','var')
            plt_val = 'b';
        else
            plt_val = plt_val.^0.4;
        end
        
        figure;
        hold on
        scatter(YR_E,XR_N,[],'^k');
        scatter(YO_E,XO_N,[],'*b');
        ddy = cos(AZI_R-pi/2)*max(dRng(:));
        ddx = -sin(AZI_R-pi/2)*max(dRng(:));
        plot([YR_E,YR_E+ddy],[XR_N,XR_N+ddx],'-.k');
        
        if scale_for_test
            scl = 4;
            azi_bin_l = azi_bin_l * scl;
            azi_bin_r = azi_bin_r * scl;
            dAzi = dAzi * scl;
        end
            
        figure;
        ax = polaraxes;
        hold on
        if size(azi_bin_l,2)==size(rng_bin_d,2)
            polarscatter(azi_bin_l,rng_bin_d,[],plt_val,'o','filled'); 
            polarscatter(azi_bin_r,rng_bin_u,[],plt_val,'x','filled');
        else
            bin_l = repmat(azi_bin_l(:),1,length(rng_bin_d(:)));
            bin_d = repmat(rng_bin_d(:),1,length(azi_bin_l(:)))';
            bin_r = repmat(azi_bin_r(:),1,length(rng_bin_u(:)));
            bin_u = repmat(rng_bin_u(:),1,length(azi_bin_r(:)))';
            polarscatter(bin_l(:),bin_d(:),[],'o');
            polarscatter(bin_r(:),bin_u(:),[],'x');
        end
        polarscatter(dAzi,dRng,[],'^k');
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir = 'clockwise';
        ax.ThetaLim = [-90 90];
        fprintf('');
        
        if scale_for_test
            azi_bin_l = azi_bin_l / scl;
            azi_bin_r = azi_bin_r / scl;
            dAzi = dAzi / scl;
        end
        
    end
    
    idx_azi = cell(1,N_points);
    idx_rng = cell(1,N_points);
    idx_azi_rng = zeros(1,N_points);
    
    for p_i = 1:N_points
        if plt_fig
            plt_ax = polarscatter(dAzi(p_i),dRng(p_i),'r','filled');
        end
        idx_azi{p_i} = find(azi_bin_l>=dAzi(p_i) & ...
                          azi_bin_r<dAzi(p_i));
        idx_rng{p_i} = find(rng_bin_d<=dRng(p_i) & ...
                          rng_bin_u>dRng(p_i));
        if isempty(idx_azi{p_i}) || isempty(idx_rng{p_i}) 
            idx_azi_rng(p_i) = 0;
        else
            if size(azi_bin_l,2) == 1
                idx_tmp = find(azi_bin_l>=dAzi(p_i) & ...
                              azi_bin_r<dAzi(p_i) & ...
                              rng_bin_d<=dRng(p_i) & ...
                              rng_bin_u>dRng(p_i));
                if ~isempty(idx_tmp)
                    idx_azi_rng(p_i) = idx_tmp(1);
                else
                    idx_azi_rng(p_i) = 0;
                end
            else
                idx_azi_rng(p_i) = ( N_theta * (idx_rng{p_i}-1) + ...
                                          idx_azi{p_i} );      
            end
        end
    end
    
    if plt_fig
        if size(azi_bin_l,2)==size(rng_bin_d,2)
            azi_bin_c = (azi_bin_l + azi_bin_r) / 2;
            rng_bin_c = (rng_bin_d + rng_bin_u) / 2;
            polarscatter(azi_bin_c(idx_azi_rng(idx_azi_rng>0)),...
                         rng_bin_c(idx_azi_rng(idx_azi_rng>0)),[],'vk');
        end
    end
    
end


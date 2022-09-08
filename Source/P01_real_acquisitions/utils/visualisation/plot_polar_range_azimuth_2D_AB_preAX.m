function [s] = plot_polar_range_azimuth_2D_AB_preAX(y_axis,...
                                              x_axis,...
                                              data,...
                                              xlimit,...
                                              ylimit,...
                                              plotfunc,...
                                              max_dev)
%PLOT_POLAR_RANGE_AZIMUTH_2D_AB Summary of this function goes here
%   Detailed explanation goes here

    if ~exist('plotfunc','var')
        plotfunc='contourf';
    end
    
    if strcmp(plotfunc,'pcolor')
        s = pcolor(y_axis,x_axis,data);
        set(s, 'edgecolor','none') 
    elseif strcmp(plotfunc,'contourf')
        if exist('max_dev','var')
            data(abs(data)>=max_dev)=NaN;
        end
        [C,s] = contourf(y_axis, x_axis, data);
        set(s,'LineColor','none');
    elseif strcmp(plotfunc,'surf')
        s = surf(y_axis, x_axis, data, 'EdgeColor','none');
    elseif strcmp(plotfunc,'scatter')
        s = scatter(y_axis(:), x_axis(:), 5, data(:),'filled');
    end
    
    view(2);
    xlabel('meters');
    ylabel('meters');
    axis equal
    if exist('xlimit','var')
        xlim(xlimit);
    else
        xlim([min(y_axis(:)),max(y_axis(:))]);
    end
    if exist('ylimit','var')
        ylim(ylimit);
    else
        ylim([min(x_axis(:)),max(x_axis(:))]);
    end
end


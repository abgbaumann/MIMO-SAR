function [keep_aoi] = get_idx_in_aoi(y_axis,x_axis, aoi)

    num_aoi = length(aoi);
    filt_aoi = zeros([length(x_axis),num_aoi]);
    
    for aoi_i = 1:num_aoi
        %% Get Bins in Area of Interest
        x_aoi = aoi{aoi_i}(:,2);
        y_aoi = aoi{aoi_i}(:,1);
        x_axis_aoi = x_axis;
        y_axis_aoi = y_axis;
        filt_aoi(:,aoi_i) = inpolygon(x_axis_aoi,y_axis_aoi,x_aoi,y_aoi);    
    end
    filter_data = sum(filt_aoi,2);
    filter_data = filter_data>0;
    keep_aoi = int32(filter_data);
    [~,keep_aoi(keep_aoi==1)] =find(filt_aoi==1); % AoI Class
    keep_aoi = logical(keep_aoi);

end


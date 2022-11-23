function [path2ts] = cascade_MIMO_03_psi2timeseries(path2proj,time_select,filt_by_asi,create_aoi)

    debugging = 0;

    path2psi = fullfile(path2proj,'03_PSI_Interfero','02_PSI.mat');
    path2ts = fullfile(path2proj,'04_TimeSeries');

    N_Hz_plt = 10;
    
    if ~isfolder(path2ts)
        mkdir(path2ts);
    end
    load(path2psi);
    

    %% Filtering by Statistic
    fprintf('data.asi\n');
    idx_asi = data.asi >= filt_by_asi{2};
    field_names = fields(data);
    for f_i = 1:length(field_names)
        if size(data.(field_names{f_i}),1) == length(idx_asi)
            data.(field_names{f_i})(~idx_asi,:) = [];
        end
    end

    %% Load Area of Interest
    path2aoi = fullfile(path2ts,'aoi.mat');
    if ~isfile(path2aoi)
        create_circles = ones(create_aoi,1) * 0;
    
        [aoi, desc] = set_aoi_scXobs(create_aoi,...
                                     data.ampl,...
                                     data.coh,...
                                     round(mean(data.coh_class,2)),...
                                     data.y_axis,...
                                     data.x_axis,...
                                     create_circles);
        save(path2aoi,'aoi','desc');
    else
        load(path2aoi);
        create_aoi = length(aoi);
    end
    
    fig = figure('units','centimeters','position',[0,0,40,20]);
    
    if time_select{1}
        idx_time = data.time_abs_interf>=time_select{2}&...
                   data.time_abs_interf<=time_select{3};
    else
        idx_time = logical(ones(size(data.time_abs_interf)));
    end

    %% Rotation
    alpha = deg2rad(90-data.radar.orientation(2)); % Elevation / Pitch
    beta = -deg2rad(data.radar.orientation(3)); % Rotation / Roll
    gamma = deg2rad(data.radar.orientation(1)); % Azimuth /Yaw
    Rx = [1,0,0;...
          0, cos(alpha), sin(alpha);...
          0, -sin(alpha), cos(alpha);];


    Ry = [cos(beta) 0 sin(beta);...
          0, 1, 0;...
          -sin(beta), 0, cos(beta)];

    Rz = [cos(gamma), -sin(gamma), 0;...
          sin(gamma), cos(gamma), 0;...
          0, 0, 1];

    CC = [data.y_axis,data.x_axis,zeros(size(data.x_axis))];

    CC_r = CC * Ry * Rx * Rz;

    data_ts.dX_North = CC_r(:,2);
    data_ts.dY_East = CC_r(:,1);
    data_ts.dZ_Height = CC_r(:,3);

    data_ts.X_North = data_ts.dX_North + data.radar.position(2);
    data_ts.Y_East = data_ts.dY_East + data.radar.position(1);
    data_ts.Z_Height = data_ts.dZ_Height + data.radar.position(3);

    data_ts.los_displ = data.cumdispl;

    data_ts.time_abs_interf = data.time_abs_interf;

    dxyz = sqrt(data_ts.dX_North.^2 + ...
                data_ts.dY_East.^2 + ...
                data_ts.dZ_Height.^2);

    %dxy = sqrt(data_ts.dX_North.^2 + data_ts.dY_East.^2);
    
    scale_displ = 1;
    data_ts.dX_North_displ = -data_ts.dX_North ./ dxyz .* data_ts.los_displ * scale_displ;
    data_ts.dY_East_displ = -data_ts.dY_East ./ dxyz .* data_ts.los_displ * scale_displ;
    data_ts.dZ_Height_displ = -data_ts.dZ_Height ./ dxyz .* data_ts.los_displ * scale_displ;

    if debugging
        fig2 = figure('units','centimeters','position',[0,0,25,25]);
        hold on
        plot3([0,1],[0,0],[0,0],'r','LineWidth',2);
        plot3([0,0],[0,1],[0,0],'r','LineWidth',2);
        plot3([0,0],[0,0],[0,1],'r','LineWidth',2);
        scatter3(0,0,0,'ok','filled');
        axis equal
        view(-8,25);
        box on
        grid on
        scatter3(CC(:,1),CC(:,2),CC(:,3),[],data.coh,'filled');
        scatter3(CC_r(:,1),CC_r(:,2),CC_r(:,3),[],data.coh,'s','filled');
    end

    if debugging
        fig3 = figure('units','centimeters','position',[0,0,25,25]);
        hold on
        scatter3(data.radar.position(1),data.radar.position(2),data.radar.position(3),'ok','filled');
        dEast = [zeros(size(data_ts.dY_East_displ,1)),data_ts.dY_East_displ(:,100)] + data_ts.Y_East;
        dNorth = [zeros(size(data_ts.dX_North_displ,1)),data_ts.dX_North_displ(:,100)] + data_ts.X_North;
        dHeight = [zeros(size(data_ts.dZ_Height_displ,1)),data_ts.dZ_Height_displ(:,100)] + data_ts.Z_Height;
        for plt_i = 1:10:size(dEast,1)
            plot3(dEast(plt_i,:),...
                     dNorth(plt_i,:),...
                     dHeight(plt_i,:),...
                     'r');
        end
        axis equal
        view(-8,25);
        box on
        grid on

        scatter3(data_ts.Y_East,data_ts.X_North,data_ts.Z_Height,[],data.coh,'filled');
    end

    figure(fig);

    data_visX = data_ts.dX_North_displ(:,idx_time);
    data_visY = data_ts.dY_East_displ(:,idx_time);
    data_visZ = data_ts.dZ_Height_displ(:,idx_time);

    data_visX = data_visX - mean(data_visX(:,1:10),2);
    data_visY = data_visY - mean(data_visY(:,1:10),2);
    data_visZ = data_visZ - mean(data_visZ(:,1:10),2);

    data_time = data.time_abs_interf(idx_time);
    
    ylimits_val = [inf, -inf];
    N_col = 6;
    sp_mat = reshape([1:N_col*create_aoi],N_col,[])';
    sp_aoi = 2;
    
    theta = data.radar.orientation(2)*pi/180;%data.radar.orientation(2)*pi/180;
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    coord = R * [data.y_axis,data.x_axis]';
    y_rot = coord(1,:)';
    x_rot = coord(2,:)';
    
    %% Overview Interferogram
    sp_id = sp_mat(:,1:sp_aoi-1);
    ax(1) = subplot(create_aoi,N_col,sp_id(:));
    hold on 
    data_map = log10(mean(data.ampl,2));
    ylimits_map = [min(y_rot),max(y_rot)];
    ylimits_map = ylimits_map + diff(ylimits_map)*0.1*[-1 1];
    xlimits_map = [min(x_rot),max(x_rot)];
    xlimits_map = xlimits_map + diff(xlimits_map)*0.1*[-1 1];
    plot_polar_range_azimuth_2D_AB_preAX(y_rot,x_rot,data_map,ylimits_map,xlimits_map,'scatter');
    box on
    grid on
    title('Amplitude [log10(A)]');
    
    for aoi_i = 1:create_aoi
        sp_id = sp_mat(aoi_i,sp_aoi:end);
        ax(aoi_i+1) = subplot(create_aoi,N_col,sp_id(:));
        hold on
        idx{aoi_i,1} = get_idx_in_aoi(data.y_axis,...
                             data.x_axis,...
                             {aoi{aoi_i}});
        data_aoiX = data_visX(idx{aoi_i},:);
        [N_bin,N_time] = size(data_aoiX);
        data_aoiX_mean = mean(data_aoiX,1);

        data_aoiY = data_visY(idx{aoi_i},:);
        data_aoiY_mean = mean(data_aoiY,1);

        data_aoiZ = data_visZ(idx{aoi_i},:);
        data_aoiZ_mean = mean(data_aoiZ,1);

        N_Hz_is = 1/seconds(data.time_rel_interf(2)-data.time_rel_interf(1));
        N_skip = round(N_Hz_is/N_Hz_plt);
        for bin_i = 1:N_bin
            plot(data_time(1:N_skip:end),data_aoiX(bin_i,1:N_skip:end),'HandleVisibility','off');
        end
        plot(data_time(1:N_skip:end),data_aoiX_mean(:,1:N_skip:end),'b','LineWidth',3,'DisplayName','Displacement to North');
        plot(data_time(1:N_skip:end),data_aoiY_mean(:,1:N_skip:end),'r','LineWidth',3,'DisplayName','Displacement to East');
        %plot(data_time(1:N_skip:end),data_aoiZ_mean(:,1:N_skip:end),'g','LineWidth',3,'DisplayName','Displacement to Height');

        ylimits_val = [min([min(ylim),ylimits_val(1)]),max([max(ylim),ylimits_val(2)])]; 
        title(desc{aoi_i});
        if aoi_i == create_aoi
            legend('Location','best');
        end

        axes(ax(1));
        coord2 = R * [aoi{aoi_i}(:,1),aoi{aoi_i}(:,2)]';
        aoi_rot{aoi_i,1} = coord2';
        pgon = polyshape(coord2(1,:),coord2(2,:));
        plot(pgon,'FaceColor','none','EdgeColor','y');

    end
    
    set(ax(1),'Color','k')
    
    for aoi_i = 1:create_aoi
        set(ax(aoi_i+1),...
            'YLim',ylimits_val,...
            'Box','on',...
            'XLim',[min(data_time),max(data_time)]);
         yline(ax(aoi_i+1),0,...
               'LineWidth',1,...
               'LineStyle','--',...
               'HandleVisibility','off');
    end
    
    Link = linkprop(ax(2:end),{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    sgtitle('Projected Displacements [mm]');
    
    path2fig = fullfile(path2ts,'Timeseries.png');
    exportgraphics(fig,path2fig,'Resolution',600);
    savefig(fig,replace(path2fig,'.png','.fig'));

    data_aoi.aoi_idx = idx;
    data_aoi.x_axis = data.x_axis;
    data_aoi.y_axis = data.y_axis;
    data_aoi.desc = desc;
    data_aoi.aoi = aoi;
    data_aoi.aoi_rot = aoi_rot;
    data_ts.coh = data.coh;
    data_ts.coh_std = data.coh_std;
    data_ts.asi = data.asi;
    data_ts.ampl = data.ampl;
    
    save(replace(path2fig,'.png','.mat'),'data_ts','data_aoi');
end


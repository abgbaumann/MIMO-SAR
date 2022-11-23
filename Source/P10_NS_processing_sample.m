close all
clear all
clc

%% Add path
[file_dir,~,~] = fileparts(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath(fullfile(file_dir)));

%% Processing Settings
plt_fig = 0; % Slow if activated.

% setting_processing.type      = {'random',...
%                                 'tilting',...
%                                 'bridge_bending',...
%                                 'sinusodial',...
%                                 'rotating',...
%                                 'bending',...
%                                 'stepping',...
%                                 }; % Processing with given Pattern
setting_processing.type      = {'tilting'};  

setting_processing.trg_noise = [2.5E-5,...
                                2.5E-5,...
                                2.5E-5,...
                                2.5E-5,...
                                2.5E-5,...
                                2.5E-5,...
                                2.5E-5,...
                                ]; % in [m], add Simulated Noise of 25 um

setting_processing.trg_disp  = [5E-3,...
                                5E-3,...
                                5E-3,...
                                5E-3,... 
                                5E-3,... 
                                5E-3,... 
                                2.0E-3,...
                                ]; % in [m], Anticipated maximum displacement
                            
setting_processing.rotations = {0,0,0,'deg'}; % Azimuth, Elevation, Rotation [degree]

scale = 1/1000; % scale adjusting for mm or m input (for Final Plotting)

scale_filt = 3; % scale_filt = 1 is equal to the TIDEP-01012 elevation Field-of-View. scale_filt = 3 = 3x initial FoV

new_random_instr = 0; % Keep predefinded Instrument location (=0) or define new random locations (>0)

N_nodisp = 10; % Number of simulated acquisitions with no displacement at the beginning and end of the movement.

N_time_sim = 25; % Number of simulated acquistions with displacement with a pattern as definded in setting_processing.type

N_stdv = N_nodisp; % Calculating the Observation Standard Deviation

initialVars = who;
initialVars{end+1} = 'initialVars';
initialVars{end+1} = 's_i';

for s_i = 1:length(setting_processing.type)

    clearvars('-except',initialVars{:})
    
    trg_noise = setting_processing.trg_noise(s_i);
    proc_type = setting_processing.type{s_i};
    trg_disp = setting_processing.trg_disp(s_i);

    %% Get Radar Positions [X, Y, Z]
    if new_random_instr == 0
        XR_N = [-0.9636    0.1632    2.0548    2.4027    2.3826];
        YR_E = [-0.3083    1.4841   -0.9021    1.2147   -0.5824];
        ZR_H = [-0.1083   -1.4781   -0.6448   -2.3553   -2.9119];
        N_radar = length(XR_N);
        AZR = zeros(1,N_radar);
    else
        N_radar = new_random_instr;
        XX = [-1,3];
        XR_N = randi(XX,1,N_radar) +  randn(1,N_radar) * 0.5;
        YY = [-4,4];
        YR_E = randi(YY,1,N_radar) + randn(1,N_radar) * 0.5;
        ZZ = [-3,0];
        ZR_H = randi(ZZ,1,N_radar) + randn(1,N_radar) * 0.5;
        AZR = zeros(1,N_radar);
    end
    r_noise = ones(1,N_radar) * trg_noise; % [m]

    for r_i = 1:N_radar
        [X_N{r_i}, Y_E{r_i}, Z_H{r_i}, ...
         x_along{r_i}, y_cross{r_i}, ...
         theta_c{r_i}, theta_l{r_i}, theta_r{r_i},...
         range_c{r_i}, range_u{r_i}, range_d{r_i}] = ...
            get_radar(XR_N(r_i),YR_E(r_i),ZR_H(r_i),AZR(r_i));
    end

    %% Get Moving Object over Time
    N_ptxy=30;
    dyz = 9;
    E_C = [-dyz/2,dyz/2,N_ptxy];
    up_min = 2;
    H_C = [0+ceil(max(ZR_H))+up_min,dyz+ceil(max(ZR_H))+up_min,N_ptxy];
    N_C = [15.5-1,N_time_sim];
    
    if strcmp(proc_type,'random')

        [E_0,H_0,N_0,de,dh,dn] = random_vector_field_deformations(0,1,E_C(3),...
                                                          0,1,H_C(3),...
                                                          0,1,N_C(2));

        dn_scale = dn ./ max(abs(dn(:))) .* trg_disp;
        de_scale = de ./ max(abs(de(:))) .* trg_disp;
        dh_scale = dh ./ max(abs(dh(:))) .* trg_disp;

        e_vector=linspace(E_C(1),E_C(2),E_C(3));
        h_vector=linspace(H_C(1),H_C(2),H_C(3));
        [ee,hh]=meshgrid(e_vector',h_vector');
        Y_0 = ee(:)';
        X_0 = repmat(N_C(1),size(Y_0));
        Z_0 = hh(:)';

        X_obj = repmat(X_0,size(dn,1),1) + dn_scale;
        Y_obj = repmat(Y_0,size(de,1),1) + de_scale;
        Z_obj = repmat(Z_0,size(dh,1),1) + dh_scale;

    elseif strcmp(proc_type,'bending')
        f = 1;
        y0 = E_C(1);
        y1 = E_C(2);
        dy = (E_C(2)-E_C(1))/E_C(3);
        X0 = N_C(1);
        N  = N_C(2);
        z0 = H_C(1);
        z1 = H_C(2);
        dz = (H_C(2)-H_C(1))/H_C(3);
        Z0 = H_C(1);

        x = @(a,y,z) a * sin(pi*f*(y-y0) / (y1-y0) ) + a * sin(pi*f*(z-z0) / (z1-z0) );
        [X_obj,Y_obj,Z_obj] = param_3d_object(y0,y1,dy,x,z0,z1,dz,-trg_disp,N,X0,plt_fig);
    elseif strcmp(proc_type,'tilting')
        y0 = E_C(1);
        y1 = E_C(2);
        dy = (E_C(2)-E_C(1))/E_C(3);
        X0 = N_C(1);
        N  = N_C(2);
        z0 = H_C(1);
        z1 = H_C(2);
        dz = (H_C(2)-H_C(1))/H_C(3);
        Z0 = H_C(1);
        x = @(a, z) a * z;
        [X_obj,Y_obj,Z_obj] = planar_3d_object(y0,y1,dy,x,z0,z1,dz,-trg_disp,N,X0,plt_fig);
    elseif strcmp(proc_type,'rotating')
        y0 = E_C(1);
        y1 = E_C(2);
        dy = (E_C(2)-E_C(1))/E_C(3);
        X0 = N_C(1);
        N  = N_C(2);
        z0 = H_C(1);
        z1 = H_C(2);
        dz = (H_C(2)-H_C(1))/H_C(3);
        Z0 = H_C(1);
        x = @(a, y) a * y;
        [X_obj,Y_obj,Z_obj] = planar_3d_object(y0,y1,dy,x,z0,z1,dz,-trg_disp,N,X0,plt_fig,'Y');
    elseif strcmp(proc_type,'stepping')
        y0 = E_C(1);
        y1 = E_C(2);
        dy = (E_C(2)-E_C(1))/E_C(3);
        X0 = N_C(1);
        N  = floor(N_C(2)/2)+N_nodisp/2;
        N_nodisp = floor(N_nodisp/2);
        z0 = H_C(1);
        z1 = H_C(2);
        dz = (H_C(2)-H_C(1))/H_C(3);
        Z0 = H_C(1);
        dx = trg_disp/5;
        [X_obj,Y_obj,Z_obj] = stepping_3d_object(y0,y1,dy,dx,z0,z1,dz,-trg_disp,N,X0,plt_fig);
        X_obj = [X_obj;flipud(X_obj)];
        Y_obj = [Y_obj;Y_obj];
        Z_obj = [Z_obj;Z_obj];
    elseif strcmp(proc_type,'sinusodial')
        y0 = E_C(1);
        y1 = E_C(2);
        dy = (E_C(2)-E_C(1))/E_C(3);
        z0 = H_C(1);
        z1 = H_C(2);
        dz = (H_C(2)-H_C(1))/H_C(3);
        Z0 = H_C(1);
        X0 = N_C(1);
        N  = N_C(2); % Number of Observations
        Lambda = 0.5; % Damping Factor
        N_waves = 5; % Number of Waves
        [X_obj,Y_obj,Z_obj] = sinusodial_time_3d_object(y0,y1,dy,z0,z1,dz,-trg_disp,N,X0,Lambda,N_waves,plt_fig);
    elseif strcmp(proc_type,'bridge_bending')
        y0 = E_C(1);
        y1 = E_C(2);
        dy = (E_C(2)-E_C(1))/E_C(3);
        z0 = H_C(1);
        z1 = H_C(2);
        dz = (H_C(2)-H_C(1))/H_C(3);
        Z0 = H_C(1);
        X0 = N_C(1);
        N  = N_C(2); % Number of Observations
        N_waves = 0.5; % Number of Waves
        [X_obj,Y_obj,Z_obj] = gaussian_spatial_3d_object(y0,y1,dy,z0,z1,dz,-trg_disp,N,X0,plt_fig);
    end
    
    %% Rotation
    alpha = setting_processing.rotations{1};
    beta =  setting_processing.rotations{2};
    gamma = setting_processing.rotations{3};
    unit =  setting_processing.rotations{4};
    [mc,nc] = size(X_obj);
    XYZ_obj_avg = [mean(X_obj(:)),mean(Y_obj(:)),mean(Z_obj(:))];
    [X_obj_tmp,Y_obj_tmp,Z_obj_tmp] = cart2cart_proj(X_obj(:)',...
                                                     Y_obj(:)',...
                                                     Z_obj(:)',...
                                                     XYZ_obj_avg(1),...
                                                     XYZ_obj_avg(2),...
                                                     XYZ_obj_avg(3),...
                                                     alpha,...
                                                     beta,...
                                                     gamma,...
                                                     unit);
    
    X_obj = reshape(X_obj_tmp,mc,nc) + XYZ_obj_avg(1);
    Y_obj = reshape(Y_obj_tmp,mc,nc) + XYZ_obj_avg(2);
    Z_obj = reshape(Z_obj_tmp,mc,nc) + XYZ_obj_avg(3);
    
    %% Outname
    str_time = datestr(datetime,'yyyymmdd_hhMMss');
    name_out = sprintf('%s_%02drad_disp%05dum_sig%03dum_a%03db%03dg%03d%s_%s',...
                            proc_type,...
                            N_radar,...
                            floor(trg_disp*1E6),...
                            floor(trg_noise*1E6),...
                            floor(alpha),...
                            floor(beta),...
                            floor(gamma),...
                            unit,...
                            str_time);

    path_out = fullfile('D00_sample_data','simulated',name_out);
    mkdir(path_out);

    %% Add Time of no deformation
    X_obj = cat(1, repmat(X_obj(1,:),N_nodisp,1),X_obj,repmat(X_obj(end,:),N_nodisp,1));
    Y_obj = cat(1, repmat(Y_obj(1,:),N_nodisp,1),Y_obj,repmat(Y_obj(end,:),N_nodisp,1));
    Z_obj = cat(1, repmat(Z_obj(1,:),N_nodisp,1),Z_obj,repmat(Z_obj(end,:),N_nodisp,1));

    %% Deformation Video
    if plt_fig
        xlim_val = [min(X_obj(:)),max(X_obj(:))];
        ylim_val = [min(Y_obj(:)),max(Y_obj(:))];
        zlim_val = [min(Z_obj(:)),max(Z_obj(:))];

        writerObj = VideoWriter(fullfile(path_out,'DefVideo_PC'));
        writerObj.Quality = 100;
        writerObj.FrameRate = 5;
        open(writerObj);

        figure('units','centimeters','position',[0,0,40,20]);
        for t_i = 1:size(Y_obj,1)
            scatter3(Y_obj(t_i,:),X_obj(t_i,:),Z_obj(t_i,:),'k','filled');
            xlim(ylim_val);ylim(xlim_val);zlim(zlim_val);
            set(gcf, 'color', [1 1 1]);
            title(sprintf('% 5d/%d',t_i,size(Y_obj,1)));
            %pause(0.2);
            writeVideo(writerObj,getframe(gcf));
        end
        axis equal
        close(gcf)
        close(writerObj);
    end
    
    %% Radar and Point Cloud
    if plt_fig
        figure('units','centimeters','position',[0,0,40,20]);
        hold on; axis equal; box on; grid on
        for r_i = 1:N_radar
            scatter3(Y_E{r_i},X_N{r_i},Z_H{r_i});
        end

        figure; 
        hold on; axis equal; box on; grid on; view(-20,30)
        scatter3(YR_E,XR_N,ZR_H,'^k','filled');
        scatter3(Y_obj(:),X_obj(:),Z_obj(:),[],'xk');
        xlabel('Y_{East}');
        ylabel('X_{North}');
        zlabel('Z_{Height}');

    end

    %% Get Bin for Points
    N_time = size(X_obj,1);
    N_points = size(X_obj,2);
    N_bins = size(X_N{1},2);

    idx_theta = cell(1,N_radar);
    idx_range = cell(1,N_radar);
    idx_tr = cell(1,N_radar);
    for r_i = 1:N_radar
        [idx_theta{r_i}, idx_range{r_i}, idx_tr{r_i}] = ...
                                            cart2bin(XR_N(r_i),...
                                                     YR_E(r_i),...
                                                     ZR_H(r_i),...
                                                     mean(X_obj,1),...
                                                     mean(Y_obj,1),...
                                                     mean(Z_obj,1),...
                                                     range_d{r_i},...
                                                     range_u{r_i},...
                                                     theta_l{r_i},...
                                                     theta_r{r_i},...
                                                     AZR(r_i)*pi/180);
    end

    if plt_fig
        for r_i=1:N_radar
            figure('units','centimeters','position',[0,0,40,20]);
            plot(Y_E{r_i},X_N{r_i},'o'); hold on
            plot(Y_E{r_i}(idx_tr{r_i}),X_N{r_i}(idx_tr{r_i}),'-xk');
            plot(Y_obj(1,:),X_obj(1,:),'-ok');
            axis equal
        end
    end

    %% Get LOS-Displacements in [mm]
    dist = zeros(N_radar,N_bins,N_time);
    disp_nF = zeros(N_radar,N_bins,N_time); % Noise Free
    list = [];
    for r_i = 1:N_radar
        ii = 1;
        for p_i = 1:N_points
            dX = X_obj(:,p_i)-XR_N(r_i);
            dY = Y_obj(:,p_i)-YR_E(r_i);
            dZ = Z_obj(:,p_i)-ZR_H(r_i);
            dXYZ = sqrt(dX.^2 + dY.^2 + dZ.^2);
            obs_i = idx_tr{r_i}(p_i);%min_idx(rad_i,obj_i);
            if obs_i > 0
                list(r_i,ii) = obs_i; ii=ii+1;
                dist(r_i,obs_i,:) = dXYZ;
                diffdXYZ = dXYZ(1)-dXYZ;
                disp_nF(r_i,obs_i,:) = diffdXYZ;
            end
        end
    end

    noise_0 = randn(size(disp_nF));
    noise_0 = repmat(r_noise',1,size(noise_0,2)) .* (noise_0 ./ std(noise_0,[],3));
    disp = disp_nF + noise_0;

    N_time = size(disp,3);

    %% Radar Simulation
    for r_i = 1:N_radar
        data_r{r_i}.x_axis = x_along{r_i};
        data_r{r_i}.y_axis = y_cross{r_i};
        data_r{r_i}.cumdispl_1N = squeeze(disp(r_i,:,:));
        data_r{r_i}.cumdispl_1N_noiseFree = squeeze(disp_nF(r_i,:,:));
        data_r{r_i}.cumdispl = data_r{r_i}.cumdispl_1N;
        data_r{r_i}.radar.position = [YR_E(r_i), XR_N(r_i), ZR_H(r_i)];
        data_r{r_i}.radar.orientation = [AZR(r_i),0,0];
        data_r{r_i}.radar.orientation_unit = 'deg';
        data_r{r_i}.time_rel_interf = seconds(0:N_time-1);
        data_r{r_i}.time_abs_interf = datetime() + data_r{r_i}.time_rel_interf;
        data_r{r_i}.mag = ones(size(data_r{r_i}.cumdispl_1N));
        data_r{r_i}.coher = ones(size(data_r{r_i}.cumdispl_1N));
        data_r{r_i}.coher_class = data_r{r_i}.coher(:,1);
        data_r{r_i}.range = round(sqrt(data_r{r_i}.x_axis.^2+data_r{r_i}.y_axis.^2),3);
        data_r{r_i}.azimuth = round(atan2(data_r{r_i}.y_axis,data_r{r_i}.x_axis),6) * 180 / pi;
        data_r{r_i}.radar.azimuth.all_unit = 'deg';
        data_r{r_i}.radar.azimuth.all = unique(data_r{r_i}.azimuth);
        data_r{r_i}.radar.range.all = unique(data_r{r_i}.range);
        data_r{r_i}.radar.range.resolution = 0.04; % [m]
        data_r{r_i}.X_inst = XR_N(r_i);
        data_r{r_i}.Y_inst = YR_E(r_i);
        data_r{r_i}.Z_inst = ZR_H(r_i);
        data_r{r_i}.az_inst = data_r{r_i}.radar.orientation(1);
        data_r{r_i}.el_inst = data_r{r_i}.radar.orientation(2);
        data_r{r_i}.rot_inst = data_r{r_i}.radar.orientation(3);
        data_r{r_i}.radar.elevation.resolution = 20; % [deg]
    end

    %% Point Cloud
    data_pc.X_north = mean(X_obj,1);
    data_pc.Y_east = mean(Y_obj,1);
    data_pc.Z_height = mean(Z_obj,1);
    data_pc.color = zeros(size(data_pc.Y_east,2),3);
    data_pc.id = 1:length(data_pc.X_north);

    %% 
    step_i = 0;
    tic0 = tic;

    initialVars = [];
    initialVars = who;
    initialVars{end+1} = 'data';
    initialVars{end+1} = 'data_r';
    initialVars{end+1} = 'data_pc';
    initialVars{end+1} = 'N_radar';
    initialVars{end+1} = 'step_i';
    initialVars{end+1} = 'rad_i';
    clearvars('-except',initialVars{:})

    for r_i = 1:N_radar

        data_pc_tmp = data_pc;

        step_i = step_i+1;
        fprintf("Step % 2d (% 6.1fs): Convert Point Cloud Coordinates to Range, Azimuth and Elevation\n",step_i,round(toc(tic0),1));
        [data_pc_tmp.x_along,data_pc_tmp.y_cross,data_pc_tmp.z_height] = ...
                                     cart2cart_proj(data_pc_tmp.X_north,...
                                                    data_pc_tmp.Y_east,...
                                                    data_pc_tmp.Z_height,...
                                                    data_r{r_i}.X_inst,...
                                                    data_r{r_i}.Y_inst,...
                                                    data_r{r_i}.Z_inst,...
                                                    data_r{r_i}.az_inst,...
                                                    data_r{r_i}.el_inst,...
                                                    data_r{r_i}.rot_inst,...
                                                    data_r{r_i}.radar.orientation_unit);

        [data_pc_tmp.range, data_pc_tmp.phi, data_pc_tmp.theta] = ...
                                     cart2polar(data_pc_tmp.x_along,...
                                                data_pc_tmp.y_cross,...
                                                data_pc_tmp.z_height,...
                                                data_r{r_i}.radar.orientation_unit); 

        %% Coarse Filtering of the Point Cloud
        step_i = step_i+1;
        fprintf("Step % 2d (% 6.1fs): Coarse Filtering of the Point Cloud\n",step_i,round(toc(tic0),1));
        dphi   = [min(data_r{r_i}.azimuth(:))-7,...
                  max(data_r{r_i}.azimuth(:))+7];
        drange = [min(data_r{r_i}.range(:)) - data_r{r_i}.radar.range.resolution,...
                  max(data_r{r_i}.range(:)) + data_r{r_i}.radar.range.resolution];

        dtheta = [-data_r{r_i}.radar.elevation.resolution * scale_filt,...
                   data_r{r_i}.radar.elevation.resolution * scale_filt];

        [idx_rough] = polar2filt(data_pc_tmp.range,...
                           data_pc_tmp.phi,...
                           data_pc_tmp.theta,...
                           drange,...
                           dphi,...
                           dtheta);
        idx_rough_log = zeros(size(data_pc_tmp.X_north));
        idx_rough_log(idx_rough) = true;

        data_pc_tmp.X_north(~idx_rough_log) = [];    
        data_pc_tmp.Y_east(~idx_rough_log) = [];  
        data_pc_tmp.Z_height(~idx_rough_log) = [];  
        data_pc_tmp.color(~idx_rough_log,:) = [];  
        data_pc_tmp.x_along(~idx_rough_log) = [];  
        data_pc_tmp.y_cross(~idx_rough_log) = [];  
        data_pc_tmp.z_height(~idx_rough_log) = [];  
        data_pc_tmp.range(~idx_rough_log) = [];  
        data_pc_tmp.phi(~idx_rough_log) = [];  
        data_pc_tmp.theta(~idx_rough_log) = [];  
        data_pc_tmp.id(~idx_rough_log) = [];

        %% Radar SLC in Range (center, up, down) and Azimuth (left, right)
        step_i = step_i+1;
        fprintf("Step % 2d (% 6.1fs): Convert Radar Coordinates to Range and Azimuth\n",step_i,round(toc(tic0),1));
        [data_r{r_i}.range_val, data_r{r_i}.phi_val,...
         data_r{r_i}.range_c, data_r{r_i}.phi_c,...
         data_r{r_i}.range_d, data_r{r_i}.range_u,...
         data_r{r_i}.phi_l, data_r{r_i}.phi_r,...
         data_r{r_i}.y_cross_axis_c,data_r{r_i}.x_along_axis_c] = ...
                                    cart2binrng2_c(data_r{r_i}.y_axis,...
                                                   data_r{r_i}.x_axis,...
                                                   data_r{r_i}.radar.azimuth.all,...
                                                   data_r{r_i}.radar.range.all,...
                                                   data_r{r_i}.radar.azimuth.all_unit,...
                                                   data_r{r_i}.radar.orientation_unit);

        [data_pc_tmp.idx_theta, data_pc_tmp.idx_range, data_pc_tmp.idx_tr] = ...
                                            cart2bin(data_r{r_i}.X_inst,...
                                                     data_r{r_i}.Y_inst,...
                                                     data_r{r_i}.Z_inst,...
                                                     data_pc_tmp.X_north,...
                                                     data_pc_tmp.Y_east,...
                                                     data_pc_tmp.Z_height,...
                                                     data_r{r_i}.range_d,...
                                                     data_r{r_i}.range_u,...
                                                     data_r{r_i}.phi_l*pi/180,...
                                                     data_r{r_i}.phi_r*pi/180,...
                                                     data_r{r_i}.az_inst*pi/180);

        %% Connect Point Cloud with Range/Azimuth-Bins
        idx_log = data_pc_tmp.idx_tr>0;

        data_pc_tmp.X_north(~idx_log) = [];          
        data_pc_tmp.Y_east(~idx_log) = [];  
        data_pc_tmp.Z_height(~idx_log) = [];  
        data_pc_tmp.color(~idx_log,:) = [];  
        data_pc_tmp.x_along(~idx_log) = [];  
        data_pc_tmp.y_cross(~idx_log) = [];  
        data_pc_tmp.z_height(~idx_log) = [];  
        data_pc_tmp.range(~idx_log) = [];  
        data_pc_tmp.phi(~idx_log) = [];  
        data_pc_tmp.theta(~idx_log) = [];  
        data_pc_tmp.idx_theta(~idx_log) = []; 
        data_pc_tmp.idx_range(~idx_log) = []; 
        data_pc_tmp.idx_tr(~idx_log) = [];
        data_pc_tmp.id(~idx_log) = [];

        data{r_i}.X_north = data_pc_tmp.X_north';
        data{r_i}.Y_east = data_pc_tmp.Y_east';
        data{r_i}.Z_height = data_pc_tmp.Z_height';
        data{r_i}.color = data_pc_tmp.color;
        data{r_i}.cumdispl = data_r{r_i}.cumdispl_1N(data_pc_tmp.idx_tr,:);  
        data{r_i}.stdv = std(data{r_i}.cumdispl(:,1:N_stdv),[],2);
        data{r_i}.idx = data_pc_tmp.id;
        data{r_i}.cumdispl_nF = data_r{r_i}.cumdispl_1N_noiseFree(data_pc_tmp.idx_tr,:);
        
        %% FOV Displacements
        dX = (data_r{r_i}.X_inst - data{r_i}.X_north);
        dY = (data_r{r_i}.Y_inst - data{r_i}.Y_east);
        dZ = (data_r{r_i}.Z_inst - data{r_i}.Z_height);

        dXYZ = sqrt(dX.^2 + dY.^2 + dZ.^2);
        dX = dX ./ dXYZ;
        dY = dY ./ dXYZ;
        dZ = dZ ./ dXYZ;

        data{r_i}.dX = repmat(dX,1,size(data{r_i}.cumdispl,2)) .* data{r_i}.cumdispl;
        data{r_i}.dY = repmat(dY,1,size(data{r_i}.cumdispl,2)) .* data{r_i}.cumdispl;
        data{r_i}.dZ = repmat(dZ,1,size(data{r_i}.cumdispl,2)) .* data{r_i}.cumdispl;

        data{r_i}.dX_nF = repmat(dX,1,size(data{r_i}.cumdispl_nF,2)) .* data{r_i}.cumdispl_nF; % Noise Free
        data{r_i}.dY_nF = repmat(dY,1,size(data{r_i}.cumdispl_nF,2)) .* data{r_i}.cumdispl_nF;
        data{r_i}.dZ_nF = repmat(dZ,1,size(data{r_i}.cumdispl_nF,2)) .* data{r_i}.cumdispl_nF;
        
        Xpos = data{r_i}.dX > 0;
        Ypos = data{r_i}.dY > 0;
        Zpos = data{r_i}.dZ > 0;
        XYZpos = (Xpos+Ypos+Zpos);
        data{r_i}.dXYZFactor = rem(XYZpos,2);
        data{r_i}.dXYZFactor(data{r_i}.dXYZFactor==0) = -1;
        data{r_i}.time_rel_interf = data_r{r_i}.time_rel_interf;
        data{r_i}.time_abs_interf = data_r{r_i}.time_abs_interf;

        clearvars('-except',initialVars{:});

    end

    %% Find all common points
    step_i = step_i+1;
    fprintf("Step % 2d (% 6.1fs): Find all overlapping observations\n",step_i,round(toc(tic0),1));
    rad_perm = perms(1:N_radar);
    rad_val = unique(sort(rad_perm(:,1:2),2),'row');
    com_idx = [];
    for comb_i = 1:size(rad_val,1)
        int_tmp = intersect(data{rad_val(comb_i,1)}.idx,data{rad_val(comb_i,2)}.idx);
        com_idx{comb_i} = int_tmp;
    end

    if size(rad_val,1)>=2
        for comb_i = 2:size(rad_val,1)
            if comb_i == 2
                idx_tmp = ismember(com_idx{comb_i-1},com_idx{comb_i});
                idx_comb = com_idx{comb_i-1}(idx_tmp);
            else
                idx_tmp = ismember(idx_comb,com_idx{comb_i});
                idx_comb = idx_comb(idx_tmp);
            end
        end
    else
        idx_comb = com_idx{1};
    end

    %% Copy and merge only common points
    step_i = step_i+1;
    fprintf("Step % 2d (% 6.1fs): Combine all overlapping observations\n",step_i,round(toc(tic0),1));
    N_points = length(idx_comb);

    data_comb.dX = cell(N_points,N_radar);
    data_comb.dY = cell(N_points,N_radar);
    data_comb.dZ = cell(N_points,N_radar);
    data_comb.dX_nF = cell(N_points,N_radar);
    data_comb.dY_nF = cell(N_points,N_radar);
    data_comb.dZ_nF = cell(N_points,N_radar);
    data_comb.dXYZFactor = cell(N_points,N_radar);


    for com_i = 1:N_points
        data_comb.X_north(com_i,1) = data{1}.X_north(data{1}.idx==idx_comb(com_i));
        data_comb.Y_east(com_i,1) = data{1}.Y_east(data{1}.idx==idx_comb(com_i));
        data_comb.Z_height(com_i,1) = data{1}.Z_height(data{1}.idx==idx_comb(com_i));

        for r_i = 1:N_radar
            if com_i == 1
                data_comb.time_rel_interf{r_i} = data{r_i}.time_rel_interf;
                data_comb.time_abs_interf{r_i} = data{r_i}.time_abs_interf;
            end
            idx_cr = data{r_i}.idx==idx_comb(com_i);
            data_comb.dX{com_i,r_i} = data{r_i}.dX(idx_cr,:);
            data_comb.dY{com_i,r_i} = data{r_i}.dY(idx_cr,:);
            data_comb.dZ{com_i,r_i} = data{r_i}.dZ(idx_cr,:);
            data_comb.dX_nF{com_i,r_i} = data{r_i}.dX_nF(idx_cr,:);
            data_comb.dY_nF{com_i,r_i} = data{r_i}.dY_nF(idx_cr,:);
            data_comb.dZ_nF{com_i,r_i} = data{r_i}.dZ_nF(idx_cr,:);
            data_comb.dXYZFactor{com_i,r_i} = data{r_i}.dXYZFactor(idx_cr,:);
            data_comb.stdv{com_i,r_i} = data{r_i}.stdv(idx_cr,:);
        end
    end

    %% Get 3D displacement Vector
    step_i = step_i+1;
    fprintf("Step % 2d (% 6.1fs): Calculate 3D displacement vector\n",step_i,round(toc(tic0),1));
    N_time = size(data_comb.time_abs_interf{1},2);
    data_comb.d_3dX = zeros(N_points,N_time);
    data_comb.d_3dY = zeros(N_points,N_time);
    data_comb.d_3dZ = zeros(N_points,N_time);
    data_comb.d_3dX_nF = zeros(N_points,N_time);
    data_comb.d_3dY_nF = zeros(N_points,N_time);
    data_comb.d_3dZ_nF = zeros(N_points,N_time);
    
    num_warn = 0;
    num_noWarn = 0;

    noise_tresh = sum(r_noise);
    warning('off');
    
    for com_i = 1:N_points
        for t_i = 1:N_time
            d_los = zeros(N_radar,3);
            d_los_nF = zeros(N_radar,3); % Noise Free
            s_los = zeros(N_radar,1);
            for r_i = 1:N_radar
                d_los(r_i,:) = [data_comb.dY{com_i,r_i}(t_i), ...
                                data_comb.dX{com_i,r_i}(t_i), ...
                                data_comb.dZ{com_i,r_i}(t_i)];
                d_los_nF(r_i,:) = [data_comb.dY_nF{com_i,r_i}(t_i), ...
                                   data_comb.dX_nF{com_i,r_i}(t_i), ...
                                   data_comb.dZ_nF{com_i,r_i}(t_i)];
                s_los(r_i,:) = data_comb.stdv{com_i,r_i};
            end
            if sum(abs(d_los(:,3))) <= noise_tresh
                d_3d_tmp = LOS23D(d_los(:,1:2), s_los);
                d_3d_tmp = [d_3d_tmp,0];
                if sum(abs(d_los_nF),'all')==0
                    d_3d_nF_tmp = [0,0,0];
                else
                    d_3d_nF_tmp = LOS23D(d_los_nF(:,1:2), 1); % [mm]
                    d_3d_nF_tmp = [d_3d_nF_tmp,0];
                end
            else
                d_3d_tmp = LOS23D(d_los, s_los); % [mm]
                if sum(abs(d_los_nF),'all')==0
                    d_3d_nF_tmp = [0,0,0];
                else
                    d_3d_nF_tmp = LOS23D(d_los_nF, 1); % [mm]
                end
            end
            [msgstr, msgid] = lastwarn;
            lastwarn('');
            if ~isempty(msgstr)
                num_warn = num_warn + 1;
            else
                num_noWarn = num_noWarn + 1;
            end
            data_comb.d_3dY(com_i,t_i) = d_3d_tmp(1);
            data_comb.d_3dX(com_i,t_i) = d_3d_tmp(2);
            data_comb.d_3dZ(com_i,t_i) = d_3d_tmp(3);
            data_comb.d_3dY_nF(com_i,t_i) = d_3d_nF_tmp(1);
            data_comb.d_3dX_nF(com_i,t_i) = d_3d_nF_tmp(2);
            data_comb.d_3dZ_nF(com_i,t_i) = d_3d_nF_tmp(3);
        end
    end
    fprintf('Number of warnings: %d/%d\n',num_warn,num_warn+num_noWarn);
    warning('on');

    data_comb.d_3dX(isnan(data_comb.d_3dX)) = 0;
    data_comb.d_3dY(isnan(data_comb.d_3dY)) = 0;
    data_comb.d_3dZ(isnan(data_comb.d_3dZ)) = 0;
    data_comb.d_3dX_nF(isnan(data_comb.d_3dX_nF)) = 0;
    data_comb.d_3dY_nF(isnan(data_comb.d_3dY_nF)) = 0;
    data_comb.d_3dZ_nF(isnan(data_comb.d_3dZ_nF)) = 0;
    
    %% Plot
    plt_LOSVec = 0;
    writerObj = VideoWriter(fullfile(path_out,'DefVideo_Rad'));
    writerObj.Quality = 100;
    writerObj.FrameRate = 5;
    open(writerObj);

    figure('units','centimeters','position',[0,0,40,20]);
    subplot(4,1,1:3)
    hold on
    min_radY = Inf;
    max_radY = -Inf;
    min_radZ = Inf;
    max_radZ = -Inf;
    for r_i=1:N_radar
        scatter3(data_r{r_i}.Y_inst,data_r{r_i}.X_inst,data_r{r_i}.Z_inst,[],'^k','filled');
        min_radY = min([min_radY,data_r{r_i}.Y_inst]);
        min_radZ = min([min_radZ,data_r{r_i}.Z_inst]);
        max_radY = max([max_radY,data_r{r_i}.Y_inst]);
        max_radZ = max([max_radZ,data_r{r_i}.Z_inst]);
    end

    scatter3(data_comb.Y_east,...
             data_comb.X_north,...
             data_comb.Z_height,...
             2,...
             'ko','filled');

    box on
    axis equal
    view(-20,25);

    dx = sqrt(data_comb.d_3dX.^2 + ...
              data_comb.d_3dY.^2 + ...
              data_comb.d_3dZ.^2);
    ddlim = [-1,1];      
    xlimits_val = [min([data_comb.Y_east(:);min_radY;YR_E(:)]),...
                   max([data_comb.Y_east(:);min_radY;YR_E(:)])] + ddlim;
    ylimits_val = ylim;
    ylimits_val = [min(ylimits_val),max(ylimits_val)+max(dx(:))/scale];
    zlimits_val = [min([data_comb.Z_height(:);min_radZ;ZR_H(:)]),...
                   max([data_comb.Z_height(:);min_radZ;ZR_H(:)])] + ddlim;

    subplot(4,1,4)
    hold on
    skip_N = 2;
    for i =1:size(dx,1)
        plot(dx(i,:),'k')
    end

    for t_i = 1:skip_N:N_time
        subplot(4,1,1:3);
        hold on
        XYZ0 = [data_comb.Y_east,...
                data_comb.X_north,...
                data_comb.Z_height];
        XYZ1 = XYZ0 + [data_comb.d_3dY(:,t_i),...
                       data_comb.d_3dX(:,t_i),...
                       data_comb.d_3dZ(:,t_i)] / scale;

        for p_i = 1:N_points
            p(p_i) = plot3([XYZ0(p_i,1),XYZ1(p_i,1)],...
                          [XYZ0(p_i,2),XYZ1(p_i,2)],...
                          [XYZ0(p_i,3),XYZ1(p_i,3)],...
                          'r-'); 
            if plt_LOSVec
                for r_i = 1:N_radar
                    dXYZ2 = [data_comb.dY{p_i,r_i}(t_i),...
                             data_comb.dX{p_i,r_i}(t_i),...
                             data_comb.dZ{p_i,r_i}(t_i)] / scale;
                    XYZ2 = XYZ0(p_i,:) + dXYZ2;
                    pp(p_i,r_i) = plot3([XYZ0(p_i,1),XYZ2(1)],...
                                  [XYZ0(p_i,2),XYZ2(2)],...
                                  [XYZ0(p_i,3),XYZ2(3)],...
                                  'b-');
                end
            end
        end

        xlim(xlimits_val);
        ylim(ylimits_val);
        zlim(zlimits_val);
        view(-45,15);
        xlabel('Y (East) [m]');
        ylabel('X (North) [m]');
        zlabel('Z (Height) [m]');
        title(sprintf('% 4d/%d',t_i,N_time))
        subplot(4,1,4);
        p(end+1) = scatter(repmat(t_i,N_points,1),dx(:,t_i),[],'r','filled');
        set(gcf, 'color', [1 1 1]);
        writeVideo(writerObj,getframe(gcf))

        delete(p);  
        if plt_LOSVec
            delete(pp); 
        end

    end

    close(gcf)
    close(writerObj);

    %% Reshape for Continues Processing

    n_t = length(data_comb.dY{1,1});
    n_pt = size(data_comb.dY,1);
    n_beob = size(data_comb.dY,2);
    for n_i = 1:n_pt
        for b_i = 1:n_beob
            XVectMat(n_i,b_i,:) = data_comb.dX{n_i,b_i};
            YVectMat(n_i,b_i,:) = data_comb.dY{n_i,b_i};
            ZVectMat(n_i,b_i,:) = data_comb.dZ{n_i,b_i};

            dXYZ_vect = [squeeze(XVectMat(n_i,b_i,:)),...
                         squeeze(YVectMat(n_i,b_i,:)),...
                         squeeze(ZVectMat(n_i,b_i,:))];
            [L(n_i,b_i,:),E(n_i,b_i,:)] = vect2eigen(dXYZ_vect);
            
            XVectMat_nF(n_i,b_i,:) = data_comb.dX_nF{n_i,b_i};
            YVectMat_nF(n_i,b_i,:) = data_comb.dY_nF{n_i,b_i};
            ZVectMat_nF(n_i,b_i,:) = data_comb.dZ_nF{n_i,b_i};
            
            dXYZ_vect_nF = [squeeze(XVectMat_nF(n_i,b_i,:)),...
                         squeeze(YVectMat_nF(n_i,b_i,:)),...
                         squeeze(ZVectMat_nF(n_i,b_i,:))];
            [L_nF(n_i,b_i,:),~] = vect2eigen(dXYZ_vect_nF);
        end
        S(n_i,:)= std(squeeze(L(n_i,:,1:N_stdv)),[],2);
        S_nF(n_i,:)= std(squeeze(L_nF(n_i,:,1:N_stdv)),[],2);
    end

    C = [data_comb.X_north, data_comb.Y_east, data_comb.Z_height];

    PC.X_obj = X_obj;
    PC.Y_obj = Y_obj;
    PC.Z_obj = Z_obj;

    path_outf = fullfile(path_out,'LECS.mat');
    save(path_outf,'L','E','C','S','L_nF','S_nF');
    path_outf = fullfile(path_out,'Comb.mat');
    save(path_outf,'data_comb');
    path_outf = fullfile(path_out,'Radar_PC.mat');
    save(path_outf,'data_pc');
    path_outf = fullfile(path_out,'Radar.mat');
    save(path_outf,'data_r');
    path_outf = fullfile(path_out,'PC_Raw.mat');
    save(path_outf,'PC');

end

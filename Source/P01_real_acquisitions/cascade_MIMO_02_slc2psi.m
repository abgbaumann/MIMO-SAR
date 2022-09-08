
%% cascade_MIMO_02_slc2cum.m
%
% The following function takes as an input the SLCs created by the function
% cascade_MIMO_01_raw2slc.m and calculates the interferogram of the time
% series.
%
% -------------------------------------------------------------------------
%
% Input:
%
% - path2proj:      A string refering to a structured folder with raw data.
%                   The directory should have a folders and a file.
%                   I) In "02_SLC_Radar_Data":
%                    a) The results created with cascade_MIMO_01_raw2slc.m
%                   II) A m-file with the same name as the proj_name and
%                   describing various acquisition properties and loading
%                   the data in the folder "02_SLC_Radar_Data".
%                   
% - filt_by_rng:    A cell consisting of three entries describing the 
%                   filtering based on slant range limitations: 
%                   {true/false, minDist, maxDist}
%                    a) The first entry defines if the distance filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the minimal distance to
%                       keep in meters.
%                    c) The third entry defines the maximal distance to
%                       keep in meters.
% - filt_by_lr:     A cell consisting of three entries describing the 
%                   filtering based on cross range limitations: 
%                   {true/false, minDist, maxDist}
%                    a) The first entry defines if the distance filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the left cross range in 
%                       meters.
%                    c) The third entry defines the right cross range in
%                       meters.
%
% - filt_by_azi:    A cell consisting of three entries describing the 
%                   filtering based on azimuth: 
%                   {true/false, minAzimuth, maxAzimuth}
%                    a) The first entry defines if the azimuth filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the minimal azimuth [deg].
%                    c) The third entry defines the maximal azimuth [deg].

% - filt_by_asi:    A cell consisting of two entries describing the 
%                   filtering based on the amplitude stability index: 
%                   {true/false, minASI}
%                    a) The first entry defines if the ASI filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the minimal ASI required.

% - filt_by_coh:    A cell consisting of two entries describing the 
%                   filtering based on the coherence [3x3 neighbourhood]: 
%                   {true/false, minCoh}
%                    a) The first entry defines if the COH filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the minimal COH required.

% - filt_by_mrd:    A cell consisting of three entries describing the 
%                   filtering based on a maximum reasonable displacement: 
%                   {true/false, maxAbsDisp, maxAvgDisp}
%                    a) The first entry defines if the MRD filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the maximal absolute 
%                       displacement [m].
%                    c) The third entry defines the maximal average
%                       displacement [m].
%
% - filt_by_time:   A cell consisting of three entries describing the 
%                   filtering based on a date/time selection: 
%                   {true/false, minDatetime, maxDatetime}
%                    a) The first entry defines if the MRD filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the start time [datetime].
%                    c) The third entry defines the end time [datetime].
%
% - filt_by_aoi:    A cell consisting of two entries describing the
%                   filtering based on N area of interests: 
%                   {true/false, numAOI}
%                    a) The first entry defines if the AOI filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the number of AoI to be
%                       expected.
%
% -------------------------------------------------------------------------
%                          Overall Workflow
%
%
%                    +------------------------------+
%   (1)              |  cascade_MIMO_01_raw2slc.m   |
%                    +--------------+---------------+
%                                   |
%                                   |
%                    +--------------¦---------------+
%   (2)              |  cascade_MIMO_02_slc2psi.m   |
%                    +--------------+---------------+  
%                                   |
%          +------------------------+------------------------+
%          |                                                 |
% 
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich (27th July 2022)


function [] = cascade_MIMO_02_slc2psi(path2proj,...
                                      filt_by_rng,...
                                      filt_by_lr,...
                                      filt_by_azi,...
                                      filt_by_asi,...
                                      filt_by_coh,...
                                      filt_by_mrd,...
                                      filt_by_time,...
                                      filt_by_aoi)

    [~,proj_name,~]= fileparts(path2proj);
    name_folder = sprintf('03_PSI_Interfero');
    path2store = fullfile(path2proj,name_folder);

    if ~isfolder(path2store)
        mkdir(path2store)
    end

    %% Data to be loaded
    step_i = 1;
    fprintf('Step % 2d: Processing Radar Data (%s)\n', step_i, proj_name);
    run(sprintf("%s.m",proj_name));
    
    scale_rad2mm = lambda / (4*pi) * 1000; % from radian to mm
    N_obs = 2*data_hz; % Average over first 2 seconds of observations and shift.
    
    %% Temporal Filtering
    if filt_by_time{1} == 1
        idx_keep = timestamp_abs >= min(filt_by_time{2:3}) & ...
                   timestamp_abs <= max(filt_by_time{2:3});
    else
        idx_keep = ones(size(timestamp_abs));
    end
    
    cplx = cplx(:,idx_keep);
    timestamp_abs = timestamp_abs(idx_keep);
    timestamp_rel = timestamp_rel(idx_keep);
    
    field_names = fields(coh);
    for field_i = 1:length(field_names)
        dims_n = ndims(coh.(field_names{field_i}));
        if dims_n == 2
            [m1,n] = size(coh.(field_names{field_i}));
            if n == length(idx_keep)
                coh.(field_names{field_i}) = coh.(field_names{field_i})(:,idx_keep);
            end
        elseif dims_n == 3
            [m1,m2,n] = size(coh.(field_names{field_i}));
            if n == length(idx_keep)
                coh.(field_names{field_i}) = coh.(field_names{field_i})(:,:,idx_keep);
            end
        end
    end
    

    %% Geometrical Filtering based on Area of Interest
    if filt_by_aoi{1}==1
        path2aoi = fullfile(path2store,'aoi.mat');
        if ~isfile(path2aoi)
            [aoi, desc] = set_aoi_scXobs(filt_by_aoi{2},...
                                         cplx,...
                                         coh.mean(:,end),...
                                         coh.class(:,end),...
                                         y_axis,...
                                         x_axis,...
                                         filt_by_aoi{3});
            save(path2aoi,"aoi","desc");
        else
            load(path2aoi,"aoi","desc");
        end
        filt_aoi = get_idx_in_aoi(y_axis,x_axis, aoi);
    else
        filt_aoi = ones(size(sct_rng));
    end

    sct_phi = round(atan2(x_axis, y_axis) - pi/2, 6);
    sct_rng = round(sqrt(x_axis.^2 + y_axis.^2), 4);
    
    if filt_by_azi{1} == 1 % Filtering based on Azimuth
        if length(filt_by_azi) == 2
            filt_azi = abs(sct_phi)*180/pi<=abs(filt_by_azi{2});
        else
            filt_azi = sct_phi*180/pi >= min(filt_by_azi{2:3}) & ...
                       sct_phi*180/pi <= max(filt_by_azi{2:3});
        end
    else
        filt_azi = ones(size(sct_rng));
    end
    
    if filt_by_rng{1} == 1 % Filtering based on Range
        filt_rng = sct_rng >= min(filt_by_rng{2:3}) & ....
                   sct_rng <= max(filt_by_rng{2:3}) ;
    else
        filt_rng = ones(size(sct_rng));
    end
    
    
    if filt_by_lr{1} == 1 % Filtering based on Cross-Range Borders
        filt_lr = y_axis >= min(filt_by_lr{2:3}) & ...
                  y_axis <= max(filt_by_lr{2:3});
    else
        filt_lr = ones(size(sct_rng));
    end
    
    
    filter_data = filt_rng & filt_azi & filt_lr & filt_aoi;
    
    cplx(~filter_data,:)=[];
    x_axis(~filter_data)=[];
    y_axis(~filter_data)=[];
    
    for field_i = 1:length(field_names)
        dims_n = ndims(coh.(field_names{field_i}));
        if dims_n == 2
            [m,n1] = size(coh.(field_names{field_i}));
            if m == length(filter_data)
                coh.(field_names{field_i})(~filter_data,:) = [];
            end
        elseif dims_n == 3
            warning('Not defined.');
        end
    end
    
    clearvars sct_rng sct_phi
    
    %% Statistical Filtering
    coh.mean = mean(coh.coh_1N,2);
    coh.std = std(coh.coh_1N,[],2);
    
    if filt_by_coh{1}
        filt_coh = coh.mean >= filt_by_coh{2};
    else
        filt_coh = ones(size(coh.mean));
    end
    
    sct_amp = abs(cplx);
    sct_std = std(sct_amp,[],2);
    sct_avg = mean(sct_amp,2);
    asi = 1-(sct_std./sct_avg); % Amplitude Stability Index
    
    if filt_by_asi{1}
        if ~exist('filt_asi','var')
            filt_asi = asi >= filt_by_asi{2};
        end
    else
        filt_asi = ones(size(cplx));
    end
    
    filter_data = filt_asi & filt_coh;
    
    cplx(~filter_data,:)=[];
    x_axis(~filter_data)=[];
    y_axis(~filter_data)=[];
    asi(~filter_data)=[];
    
    for field_i = 1:length(field_names)
        dims_n = ndims(coh.(field_names{field_i}));
        if dims_n == 2
            [m,n1] = size(coh.(field_names{field_i}));
            if m == length(filter_data)
                coh.(field_names{field_i})(~filter_data,:) = [];
            end
        elseif dims_n == 3
            warning('Not defined.');
        end
    end
    
    clearvars stc_amp sct_std sct_avg
    
    %% Number of Scatter / Observations
    [n_sct,n_obs] = size(cplx);
       
    %% Interferograms
    step_i = step_i+1;
    fprintf("Step % 2d: Compute Interferograms\n",step_i); % Averaging 
    data.interf = complex(zeros(n_sct,n_obs-1));
    cplx_avg = mean(cplx(:,1:N_obs),2);
    for obs_i = 1:n_obs-1
        data.interf(:,obs_i) = cplx(:,obs_i+1) .* conj(cplx_avg);
    end
    
    %% Store Time Stamp to data
    data.time_rel_interf = timestamp_rel(1:end-1) + ...
                                diff(timestamp_rel)/2; % Duration
    data.time_abs_interf = timestamp_abs(1:end-1) + ...
                                diff(timestamp_abs)/2; % Datetime
    
    data.time_rel_interf_unit = 'duration';
    data.time_abs_interf_unit = 'datetime';
    
    %% Saving SLC file
    fname = '01_SLC_Filt';
    path2mat = join([fullfile(path2store,fname),'.mat'],'');
    
    coher = coh.mean;
    coher_class = coh.class;
    coher_std = coh.std;
    class_id = coh.class_id;
    class_descr = coh.class_descr;
    
    step_i = step_i+1;
    fprintf("Step % 2d: Saving SLC file as *.mat\n",step_i);
    save(path2mat,...
        'cplx',...
        'coher',...
        'coher_class',...
        'coher_std',...
        'class_id',...
        'class_descr',...
        'y_axis','x_axis',...
        'timestamp_rel',...
        'timestamp_abs');
     
    clearvars cplx coher coher_class coher_std class_id class_descr
    
    %% Interferometric Phases
    step_i = step_i+1;
    fprintf("Step % 2d: Compute Interferometric Phases and Unwrap\n",step_i)
    data.ampl = abs(data.interf); % Magnitude of Complex Values
    data.phase = angle(data.interf); % Interferometric Phase
    
    %% Cumulative Displacement
    step_i = step_i+1;
    fprintf("Step % 2d: Compute Displacements in Line of Sights\n",step_i)
    data.cumdispl = unwrap(data.phase,[],2) * scale_rad2mm;
    
    %% Store Geometric Properties to data
    step_i = step_i+1;
    fprintf("Step % 2d: Range and Azimuth\n",step_i)
    
    data.x_axis = x_axis;
    data.x_axis_unit = 'm';
    data.y_axis = y_axis;
    data.y_axis_unit = 'm';
    data.azimuth = round(atan2(y_axis,x_axis)*180/pi,6);
    data.azimuth_unit = 'deg';
    data.range = round(sqrt(x_axis.^2+y_axis.^2),3);
    data.range_unit = 'm';
    data.asi = asi; % Amplitude Stability Index
    data.coh = coh.mean; % Coherence in [3x3] neighbourhood
    data.coh_class = coh.class; % Coherence Classes
    data.coh_std = coh.std; % Standard Deviation
    data.coh_class_id = coh.class_id; % Class ID [-2...1]
    data.coh_class_descr = coh.class_descr; % Description of Class ID
    
    %% Store Filtering Properties to data
    data.proc_filt.type = [];
    data.proc_filt.type_val = [];
    data.proc_filt.type_unit = [];
    
    if filt_by_rng{1} == 1
        data.proc_filt.type{end+1} = 'Filter by Slant Range [min/max]';
        data.proc_filt.type_val{end+1} = filt_by_rng{2:3};
        data.proc_filt.type_unit{end+1} = 'm';
    end
    
    if filt_by_lr{1} == 1
        data.proc_filt.type{end+1} = 'Filter by Cross Range [left/right]';
        data.proc_filt.type_val{end+1} = filt_by_lr{2:3};
        data.proc_filt.type_unit{end+1} = 'm';
    end
    
    if filt_by_azi{1} == 1
        data.proc_filt.type{end+1} = 'Filter by Azimuth [min/max]';
        data.proc_filt.type_val{end+1} = filt_by_azi{2:3};
        data.proc_filt.type_unit{end+1} = 'deg';
    end
    
    if filt_by_coh{1} == 1
        data.proc_filt.type{end+1} = 'Filter by Coherence [min]';
        data.proc_filt.type_val{end+1} = filt_by_coh{2};
        data.proc_filt.type_unit{end+1} = '-';
    end
    
    if filt_by_asi{1} == 1
        data.proc_filt.type{end+1} = 'Filter by Amplitude Stability Index [min]';
        data.proc_filt.type_val{end+1} = filt_by_asi{2};
        data.proc_filt.type_unit{end+1} = '-';
    end
    
    if filt_by_mrd{1} == 1
        data.proc_filt.type{end+1} = 'Filter by Maximum Displacement [MaxAbsDispl, MaxMeanDispl]';
        data.proc_filt.type_val{end+1} = filt_by_mrd{2:3};
        data.proc_filt.type_unit{end+1} = 'mm';
    end
    
    if filt_by_time{1} == 1
        data.proc_filt.type{end+1} = 'Filter by Time [minDate, maxDate]';
        data.proc_filt.type_val{end+1} = filt_by_time{2:3};
        data.proc_filt.type_unit{end+1} = 'datetime';
    end

    if filt_by_aoi{1} == 1
        data.proc_filt.type{end+1} = 'Filter by Area of Interest [y_axis, x_axis]';
        data.proc_filt.type_val{end+1} = aoi;
        data.proc_filt.type_unit{end+1} = 'm';
    end
    
    %% Store Radar Properties to data
    data.radar.position = radar_position; %[m]
    data.radar.position_unit = 'm'; 
    data.radar.orientation = radar_orient; %[degree]
    data.radar.orientation_unit = 'deg'; 
    data.radar.azimuth.num_virt_ant = num_antenna; %[-]
    data.radar.elevation.resolution = round((1/length_antenna)*180/pi,6); %[degree]
    data.radar.elevation.resolution_unit = 'deg'; 
    data.radar.range.resolution = rng_res; % [m]
    data.radar.range.resolution_unit = 'm';
    data.radar.azimuth.all = round(unique(data.azimuth),6);
    data.radar.azimuth.all_unit = 'deg'; 
    data.radar.range.all = unique(data.range);
    data.radar.range.all_unit = 'm'; 
    
    %% Store Processed Data
    step_i = step_i+1;
    fprintf("Step % 2d: Saving Mat-File.\n",step_i)
    
    fname = '02_PSI';
    path2mat = join([fullfile(path2store,fname),'.mat'],'');
    save(path2mat,'data');
    
    %% Plot Overview
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    
    ylimits = [min(data.y_axis),max(data.y_axis)];
    ylimits = ylimits + [-1 1] * abs(diff(ylimits)) * 0.05;
    xlimits = [min(data.x_axis),max(data.x_axis)];
    xlimits = xlimits + [-1 1] * abs(diff(ylimits)) * 0.05;
    
    %% Amplitude
    ax(1)=subplot(2,2,1); 
    ampl = log10(mean(data.ampl,2));
    plot_polar_range_azimuth_2D_AB_preAX(data.y_axis,...
                                         data.x_axis,...
                                         ampl,...
                                         ylimits,xlimits,'scatter');
    hcb = colorbar;
    hcb.Label.String = "Amplitude [log10(A)]";
    climits_val = [floor(prctile(ampl,50)),ceil(prctile(ampl,95))];
    clim(climits_val);
    caxis(climits_val);
    box on
    axis equal
    
    %% Amplitude Stability Index
    ax(2)=subplot(2,2,2);
    plot_polar_range_azimuth_2D_AB_preAX(data.y_axis,...
                                         data.x_axis,...
                                         data.asi,...
                                         ylimits,xlimits,'scatter');
    hcb = colorbar;
    hcb.Label.String = "Amplitude Stability Index [-]";
    climits_val = [floor(prctile(data.asi,50)*20)/20,1];
    clim(climits_val);
    caxis(climits_val);
    box on
    axis equal
    
    %% Coherence
    ax(3)=subplot(2,2,3);
    plot_polar_range_azimuth_2D_AB_preAX(data.y_axis,...
                                         data.x_axis,...
                                          data.coh,...
                                          ylimits,xlimits,'scatter');
    hcb = colorbar;
    hcb.Label.String = "Mean Coherence [0...1]";
    climits_val = [prctile(data.coh,20),prctile(data.coh,98)];
    clim(climits_val);
    caxis(climits_val);
    box on
    axis equal
    
    %% Maximal Displacement in Time
    ax(4)=subplot(2,2,4);
    prc_tile_coh = prctile(data.coh,95);
    max_abs_displ = max(abs(data.cumdispl(data.coh>=min([prc_tile_coh,0.95]),:)),[],1);
    idx_time_max = find(max_abs_displ==max(max_abs_displ));
    displace = squeeze(data.cumdispl(:,idx_time_max));
    
    plot_polar_range_azimuth_2D_AB_preAX(data.y_axis,...
                                         data.x_axis,...
                                         displace,...
                                         ylimits,xlimits,'scatter');
    hcb = colorbar;
    hcb.Label.String = "Displacement [mm]";
    climits = [prctile(displace,5),prctile(displace,50)];
    climits = [-max(abs(climits)),max(abs(climits))];
    clim(climits);
    caxis(climits);
    box on
    axis equal
    
    Link = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    
    sgtitle(sprintf('%s - Overview',replace(proj_name,'_',' ')));
    
    exportgraphics(fig,replace(path2mat,'.mat','.png'),'Resolution',600);
    savefig(fig,replace(path2mat,'.mat','.fig'));

end
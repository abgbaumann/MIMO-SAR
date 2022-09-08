%  Copyright (C) 2018 Texas Instruments Incorporated - http://www.ti.com/
%
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions
%   are met:
%
%     Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%
%     Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the
%     distribution.
%
%     Neither the name of Texas Instruments Incorporated nor the names of
%     its contributors may be used to endorse or promote products derived
%     from this software without specific prior written permission.
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
%   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%

%% cascade_MIMO_01_raw2slc.m
%
% The following scrip/function was modified from the above mentioned 
% TI-scirpt. It can be used to process raw data acquired with a TIDEP-01012 
% device. The result of this function is a dataset consisting of SLCs.
%
% -------------------------------------------------------------------------
%
% Input:
%
% - path2proj:      A string refering to a structured folder with raw data.
%                   The directory should have two folders.
%                   I) In "00_Parameter_Files":
%                    a) A sensor calibration file (Sensor_Calibration_at_77GHz.mat)
%                    b) A parameter file (module_param_1Chirp.m)
%                   II) In "01_Raw_Radar_Data":
%                    a) All acquisition files exported from the TIDEP-01012
%                       SSD (..._data.bin, ..._idx.bin)
%                    b) The json-Files describing the acquisition/settings
%                       with the ending .mmwave.json and .setup.json.
%                   
% - filt_by_dist:   A cell consisting of three entries describing the 
%                   filtering based on slant range: 
%                   {true/false, minDist, maxDist}
%                    a) The first entry defines if the distance filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the minimal distance to
%                       keep in meters.
%                    c) The third entry defines the maximal distance to
%                       keep in meters.
%
% - filt_by_azi:    A cell consisting of three entries describing the 
%                   filtering based on azimuth: 
%                   {true/false, minAzimuth, maxAzimuth}
%                    a) The first entry defines if the azimuth filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the minimal azimuth [deg].
%                    c) The third entry defines the maximal azimuth [deg].
%
% - filt_by_asi:    A cell consisting of two entries describing the 
%                   filtering based on the amplitude stability index: 
%                   {true/false, minASI}
%                    a) The first entry defines if the ASI filtering
%                       is applied (true) or not (false).
%                    b) The second entry defines the minimal ASI required.
%
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
% by Andreas Baumann-Ouyang, ETH Zürich (26th July 2022)

function [] = cascade_MIMO_01_raw2slc(path2proj,filt_by_dist,filt_by_azi,filt_by_asi)
    
    save_static_only = 1;     % Save only static data
    dataPlatform = 'TDA2';

    [~,proj_name,~]= fileparts(path2proj);
    
    path_parm = fullfile(path2proj,'00_Parameter_Files'); % Path to Processing Settings
    path_calib = fullfile(path_parm,'phaseMismatchCalibration.mat'); % Path to calibration file, for each board and chirp setting unique (!)
    path_param = char(fullfile(path_parm,'module_param_1Chirp.m')); % Path to parameter file for the 1 chirp configuration
    
    path_in = fullfile(path2proj,'01_Raw_Radar_Data'); % Input Path to Raw Data
    path_out = fullfile(path2proj,'02_SLC_Radar_Data'); % Output Path to SLC Data
    
    if ~isfolder(path_out)
        mkdir(path_out);
    end
    
    if ~isfolder(path_parm)
        mkdir(path_parm);
    end
    
    %% Start Processing
    step_i = 1; % Numbering for command window message
    fprintf('Step %02d: Process Radar Data (%s)\n', step_i, proj_name);
    
    %% Create List of Input Files and Store as TXT
    input_data_list = fullfile(path_parm,"Input_Data.txt");
    fidList = fopen(input_data_list,'w');
    
    path_files = char([path_in,filesep]);
    fprintf(fidList, "%s\n",path_files);
    fprintf(fidList, "%s\n",path_calib);
    fprintf(fidList, "%s\n",path_param);
    
    ID = 1;
    
    %% Parameter File
    para_name = sprintf('parameters_%02d.m',ID);
    pathGenParaFile = fullfile(path_parm,para_name);
    fid = fopen(pathGenParaFile,"w");
    clear(pathGenParaFile);
    
    %generate parameter file
    parameter_file_gen_json(path_files, ...
                            path_calib, ...
                            path_param, ...
                            pathGenParaFile, ...
                            dataPlatform);
    
    %load calibration parameters
    load(path_calib);
    
    % simTopObj is used for top level parameter parsing and data 
    % loading and saving
    simTopObj           = simTopCascade('pfile', pathGenParaFile);
    calibrationObj      = calibrationCascade(...
                                'pfile', pathGenParaFile, ...
                                'calibrationfilePath', path_calib);
    detectionObj        = CFAR_CASO('pfile', pathGenParaFile);
    
    totNumFrames        = simTopObj.totNumFrames;
    rangeBinSize        = detectionObj.rangeBinSize;
    antenna_azimuthonly = detectionObj.antenna_azimuthonly;
    frameCountGlobal    = 0;
    
    
    % Get Unique File Idxs in the "dataFolder_test"   
    [fileIdx_unique]    = getUniqueFileIdx(path_files);
    
    t0 = tic;
    
    fclose(fid);
    
    for i_file = 1:(length(fileIdx_unique))
    
        % Get File Names for the Master, Slave1, Slave2, Slave3   
        [fileNameStruct]= getBinFileNames_withIdx(path_files, fileIdx_unique{i_file});        
    
        % pass the Data File to the calibration Object
        calibrationObj.binfilePath = fileNameStruct;
    
        % Get Valid Number of Frames 
        [numValidFrames, ~] = getValidNumFrames(fullfile(path_files, fileNameStruct.masterIdxFile));
        
        
        % intentionally skip the first frame due to TDA2 
        if i_file == 1 % Skip the first frame of the dataset
            startValidFrame = 2;
        else
            startValidFrame = 1;
            
        end
    
        no_of_local_acquisitions = numValidFrames-(startValidFrame-1);
    
        step_i = step_i + 1;
    
        for frameIdx = startValidFrame:numValidFrames
    
            % read and calibrate raw ADC data            
            calibrationObj.frameIdx = frameIdx;
            frameCountGlobal = frameCountGlobal+1;
            frameCountLocal = frameIdx-(startValidFrame-1);
            adcData = datapath(calibrationObj);
    
            % RX Channel re-ordering
            adcData = adcData(:,:,calibrationObj.RxForMIMOProcess,:);            
    
            filt_type = 'hanning';
            numAngles = 256;
            [complex_data_static_tmp,...
             complex_data_dynamic_tmp,...
             x_axis,...
             y_axis] = adc2slc(adcData,...
                                 rangeBinSize,...
                                 filt_type,...
                                 numAngles,...
                                 antenna_azimuthonly);
    
            if frameCountLocal == 1
                complex_data_static = zeros([size(complex_data_static_tmp),no_of_local_acquisitions]);
                %complex_data_dynamic = zeros([size(complex_data_dynamic_tmp),no_of_local_acquisitions]);
            end
            
            complex_data_static(:,:,frameCountLocal) = complex_data_static_tmp;
    
            time_i = toc(t0)/frameCountGlobal;
            time_rem = time_i*(totNumFrames-frameCountGlobal);
            time_rem_min = floor(time_rem/60);
            time_rem_sec = floor(time_rem - (time_rem_min*60));
            time_end = datetime() + seconds(time_rem);
    
            if exist('lineLength','var')
                fprintf(repmat('\b',1,lineLength));
            end
    
            print_str = sprintf("Step %02d: Frame % 5d processed (SLC).",...
                                step_i,...
                                frameCountGlobal);
            if frameCountGlobal/totNumFrames<0.1
                lineLength = fprintf("%s\n",print_str);
            else
                lineLength = fprintf("%s (Ending in %d:%02dmin at %s)\n",...
                        print_str,...
                        time_rem_min,...
                        time_rem_sec,...
                        datestr(time_end,"dd.mm.yyyy HH:MM:SS"));
            end
        end
    
        clearvars lineLength
    
        ind = strfind(path_files, filesep);
        testName = path_files(ind(end-1)+1:(ind(end)-1));
        path_mat_03 = fullfile(path_out,sprintf("%s_%05d.mat",...
                                                    testName,...
                                                    i_file)); 
    
        %% Calculate Coherence of first to second, to N/2, and to N
        step_i = step_i + 1;
        fprintf('Step %02d: Calculate Coherences.\n',step_i);
            
        step_i = step_i + 1;
        coherence = zeros(size(complex_data_static));
        N_slc = size(complex_data_static,3);
        for slc_i = 1:N_slc
            coherence(:,:,slc_i) = get_coherence(complex_data_static(:,:,1),...
                                            complex_data_static(:,:,slc_i),...
                                            3);
    
            if exist('lineLength','var')
                fprintf(repmat('\b',1,lineLength));
            end
    
            print_str = sprintf("Step %02d: Frame % 5d / %d processed (Coherence).",...
                                step_i,...
                                slc_i,...
                                N_slc);
            lineLength = fprintf("%s\n",print_str);
        end
       
        coh_mean = mean(coherence,3);
        coh_std = std(coherence,[],3);
    
        coh_class_tresh = 0.9;
        coh_mat_e_filt = ones(size(coh_mean));
        coh_mat_e_filt(squeeze(coherence(:,:,end))<=coh_class_tresh) = 0;
        coh_mat_s_filt = 2 * ones(size(coh_mean));
        coh_mat_s_filt(squeeze(coherence(:,:,2))<=coh_class_tresh) = 0;
        coh_class = coh_mat_e_filt-coh_mat_s_filt;
    
        %% Visualise Data Amplitudes and Coherences
        step_i = step_i + 1;
        fprintf('Step %02d: Visualise Amplitudes and Coherences\n',step_i);
        close all
    
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        amp = abs(complex_data_static);
    
        % Amplitude Image
        ax(1)=subplot(2,2,1);
        plot_polar_range_azimuth_2D_AB_preAX(y_axis,...
                                                x_axis,...
                                                log10(mean(amp,3)));
        hcb = colorbar;
        hcb.Label.String = "Amplitude [log10(A)]";
        
        ax(2)=subplot(2,2,2);
        dim=ndims(complex_data_static);
        n_max = 4000;
        if size(complex_data_static,dim)>n_max
            n_skip = ceil(size(complex_data_static,dim)/n_max);
        else
            n_skip = 1;
        end
        if dim==3
            amptmp = squeeze(amp(:,:,1:n_skip:end));
        else
            error('Data format for Amplitude with given Dimension is unknown!\n');
        end
        amp_sigma = std(amptmp,[],dim);
        amp_mu = mean(amptmp,dim);
        plt_asi = 1 - (amp_sigma ./ amp_mu);
        plt_asi(plt_asi<=0) = 0;
        plot_polar_range_azimuth_2D_AB_preAX(y_axis,...
                                              x_axis,...
                                              plt_asi);
        hcb = colorbar;
        hcb.Label.String = "Amplitude Stability Index [-]";
        climits = [prctile(plt_asi(:),50),prctile(plt_asi(:),95)];
        clim(climits);
        caxis(climits);
    
        ax(3)=subplot(2,2,3);
        plot_polar_range_azimuth_2D_AB_preAX(y_axis,...
                                              x_axis,...
                                              coh_mean);
        hcb = colorbar;
        hcb.Label.String = "Mean Coherence [0...1]";
        climits = [prctile(coh_mean(:),50),prctile(coh_mean(:),90)];
        clim(climits);
        caxis(climits);
        
        ax(4)=subplot(2,2,4);
        plot_polar_range_azimuth_2D_AB_preAX(y_axis,...
                                              x_axis,...
                                              coh_std);
        colormap(ax(4), flipud(colormap))
        hcb = colorbar;
        hcb.Label.String = "Standard Deviation of Coherence";
        climits = [prctile(coh_std(:),5),prctile(coh_std(:),50)];
        clim(climits);
        caxis(climits);
    
        Link = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim'});
        setappdata(gcf, 'StoreTheLink', Link);
    
        sgtitle(sprintf('%s - Amplitude and Coherence (%d)',replace(proj_name,'_',' '),i_file));
        
        exportgraphics(fig,replace(path_mat_03,'.mat','_1.png'),'Resolution',600);
        savefig(fig,replace(path_mat_03,'.mat','_1.fig'));
    
        %% Saving Data to MAT-File.
        step_i = step_i + 1;
        fprintf('Step %02d: Compress Data by Filtering based:\n',step_i);
        
        if filt_by_dist{1}
            fprintf('- Range Limits (%.2f - %.2f m)\n',filt_by_dist{2},filt_by_dist{3});
            rng = sqrt(x_axis(1,:).^2+y_axis(1,:).^2);
            idx_rmv = rng<filt_by_dist{2} | rng>filt_by_dist{3};
            x_axis(:,idx_rmv) = []; 
            y_axis(:,idx_rmv) = []; 
            complex_data_static(:,idx_rmv,:) = [];
            coherence(:,idx_rmv,:) = [];
            coh_mean(:,idx_rmv) = [];
            coh_std(:,idx_rmv) = [];
            coh_class(:,idx_rmv) = [];
        end
    
        if filt_by_azi{1}
            fprintf('- Azimuth Limits (%.1f - %.1f deg)\n',filt_by_azi{2},filt_by_azi{3});
            azi = atan2(y_axis(:,1),x_axis(:,1)) * 180 / pi;
            idx_rmv = azi<filt_by_azi{2} | azi>filt_by_azi{3};
            x_axis(idx_rmv,:) = []; 
            y_axis(idx_rmv,:) = []; 
            complex_data_static(idx_rmv,:,:) = [];
            coherence(idx_rmv,:,:) = [];
            coh_mean(idx_rmv,:) = [];
            coh_std(idx_rmv,:) = [];
            coh_class(idx_rmv,:) = [];
        end
        
        x_axis_coh = x_axis;
        y_axis_coh = y_axis;
    
        if filt_by_asi{1} == 1
    
            if length(filt_by_asi{2}) == 1
                fprintf('- Amplitude Stability Index (>=%.2f)\n',filt_by_asi{2});
            end

             [complex_data_static,x_axis,y_axis,filt_by_asi{2},idx_a,idx_b] = ...
                filter_by_asi(complex_data_static,...
                              x_axis,...
                              y_axis,...
                              i_file,...
                              filt_by_asi{2});
                        
            if save_static_only~=1
                [complex_data_dynamic,~,~,~] = ...
                    filter_by_asi(complex_data_dynamic,...
                                  x_axis_old,...
                                  y_axis_old,...
                                  i_file,...
                                  filt_by_asi{2});
            end
        
            % Update Coherences
            coherence_tmp = zeros(size(complex_data_static));
            coh_mean_tmp = zeros(length(idx_a),1);
            coh_std_tmp = coh_mean_tmp;
            coh_class_tmp = coh_mean_tmp;
            for idx_i = 1:length(idx_a)
                coherence_tmp(idx_i,:) = coherence(idx_a(idx_i),idx_b(idx_i),:);
                coh_mean_tmp(idx_i,:) = coh_mean(idx_a(idx_i),idx_b(idx_i));
                coh_std_tmp(idx_i,:) = coh_std(idx_a(idx_i),idx_b(idx_i));
                coh_class_tmp(idx_i,:) = coh_class(idx_a(idx_i),idx_b(idx_i));
            end
    
            coherence = coherence_tmp;
            coh_mean = coh_mean_tmp;
            coh_std = coh_std_tmp;
            coh_class = coh_class_tmp;
            x_axis_coh = x_axis;
            y_axis_coh = y_axis;
        else
            [m,n,o] = size(complex_data_static);
            complex_data_static = reshape(complex_data_static,m*n,o);
            x_axis = reshape(x_axis,m*n,[]);
            y_axis = reshape(y_axis,m*n,[]);
            coherence = reshape(coherence,m*n,o);
            coh_mean = reshape(coh_mean,m*n,[]);
            coh_std = reshape(coh_std,m*n,[]);
            coh_class = reshape(coh_class,m*n,[]);
            x_axis_coh = x_axis;
            y_axis_coh = y_axis;
        end
    
    
    
        step_i = step_i + 1;
        fprintf('Step %02d: Visualise Filtered Image: Amplitude\n',step_i);
    
        fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
        amp = abs(complex_data_static);
        r=sqrt(y_axis.^2+x_axis.^2)*1.4*pi/180*100;
        scatter(y_axis,x_axis,r,log10(mean(amp,2)),'filled');
        hcb = colorbar;
        hcb.Label.String = "Amplitude [log10(A)]";
        axis equal
        box on
        xlabel('Cross Range [m]');
        ylabel('Along Range [m]');
        sgtitle(sprintf('%s - Filtered Amplitude Image (%d)',replace(proj_name,'_',' '),i_file));
    
        exportgraphics(fig2,replace(path_mat_03,'.mat','_2.png'),'Resolution',600);
        savefig(fig2,replace(path_mat_03,'.mat','_2.fig'));
         
    
        %% Saving Data to MAT-File
        coh.class = coh_class;
        coh.class_id = [-2,-1,0,1];
        coh.class_descr = {'Higher Coherence at Start',...
                                 'High Coherence',...
                                 'Low Coherence',...
                                 'Higher Coherence at End'};
    
        coh.coh_1N = coherence;
        coh.mean   = coh_mean;
        coh.std    = coh_std;
        coh.x_axis = x_axis_coh;
        coh.y_axis = y_axis_coh;
        
        step_i = step_i + 1;
        fprintf('Step %02d: Saving Data to MAT-File\n',step_i);
        if save_static_only
            save(path_mat_03,...
                'complex_data_static',...
                'coh',...
                'x_axis',...
                'y_axis',...
                'numAngles',...
                '-v7.3');
        else
            save(path_mat_03,...
                'complex_data_static',...
                'complex_data_dynamic',...
                'coh',...
                'x_axis',...
                'y_axis',...
                'numAngles',...
                '-v7.3');
        end
    
    end
    
    ID = ID + 1;

end

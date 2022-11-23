%% P03_LSQ_processing_sample.m
%
% The following script can be used to process real or simulated (MIMO-)SAR
% acquistions to derive 3D displacement vectors. The least square
% adjustment is carried out using parametric and non-parametric functions
% on the temporal or spatial domain.
%
% Workflow for Simulated Data:
%
% 1.) P02_NS_processing_sample.m
% 2.) P03_LSQ_processing_sample.m
%
%
% Workflow for Real Data:
%
% 1.) tbd
% 2.) P03_LSQ_processing_sample.m
% 
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich (26th July 2022)
%


clearvars
close all
clc

create_vis = 0; % Slow if activated.

%% Add paths
[file_dir,~,~] = fileparts(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath(fullfile(file_dir)));

%% Initalisation with declaration of path to input data
path_main = fullfile('D00_sample_data','simulated'); % Main Input Path for Simulated Data

files = dir(path_main);
dirFlags = [files.isdir];
path_struct = files(dirFlags);
path_files = {};
for fi = 1:length(path_struct) % Loads every folder
    name_i = path_struct(fi).name;
    if ~strcmp(name_i,'.') && ~strcmp(name_i,'..')
        path_files{end+1,1} = name_i;
    end
end

initialVars = who;
initialVars{end+1} = 'initialVars';
initialVars{end+1} = 'd_i';

for d_i =  1:length(path_files)
    
    clearvars('-except',initialVars{:});
    
    t_d0 = tic;
    fprintf('---------------------------------------------------------\n');
    fprintf('Started processing [%s]\n',path_files{d_i});
    path_in = fullfile(path_main,path_files{d_i});
    path_in_lecs = fullfile(path_in,'LECS.mat');
    path_in_comb = fullfile(path_in,'Comb.mat');

    if ~isfile(path_in_lecs) || ~isfile(path_in_comb) 

        warning('No MAT-File in Directory [%s]. Directory skipped.\n',path_in)

    else

        load(path_in_lecs); % Load L/E/C/S-Matrices
        load(path_in_comb); % Load processed data of simplified LSQ
    
        path_out = fullfile(path_in,'Advanced'); % Path for Saving Files
        if ~isfolder(path_out)
            mkdir(path_out);
        end
    
        %% Initial Variables      
        D_adj = {};
        L_adj = {};
        fileName = {};
        Desc = {};
        
        %% Parametric-Temporal LSQ Adjustment
        [D_adj{end+1}, L_adj{end+1}] = LSQ_Temp_Para(L,E,S,{'poly03'});
        fileName{end+1} = 'temporal_para_poly03';
        Desc{end+1} = 'Temporal cubic parametric LSQ';
        
        %% Non-parametric-Temporal LSQ Adjustment
        cov = [1,20;...
               1,20;...
               1,20]; % Parameters
           
        [D_adj{end+1}, L_adj{end+1}] = LSQ_Temp_Spline(L,E,S,cov);
        fileName{end+1} = 'temporal_spline';
        Desc{end+1} = sprintf('Temporal non-parametric LSQ (%.2f-%.2f)',cov(1,1),cov(1,2));
        
        %% Parametric-Spatial LSQ Adjustment   
        input = {'poly03',C};
        [D_adj{end+1}, L_adj{end+1}] = LSQ_Spatial_Para(L,E,S,input);
        fileName{end+1} = 'spatial_para_poly03';
        Desc{end+1} = 'Spatial cubic parametric LSQ';
        
        %% Non-parametric-Spatial LSQ Adjustment
        cov = [1,20;...
               1,20;...
               1,20];  % Parameters
    
        [D_adj{end+1}, L_adj{end+1}] = LSQ_Spatial_Spline(L,E,C,S,cov);
        fileName{end+1} = 'spatial_spline';
        Desc{end+1} =  sprintf('Spatial non-parametric LSQ (%.2f-%.2f)',cov(1,1),cov(1,2));
        

        [D_adj{end+1}, L_adj{end+1}] = LSQ_Spatial_Spline2(L,E,C,S,cov);
        fileName{end+1} = 'spatial_spline';
        Desc{end+1} =  sprintf('Spatial non-parametric LSQ2 (%.2f-%.2f)',cov(1,1),cov(1,2));

        %% Store Results
        path_mat = fullfile(path_out,...
                            sprintf('LSQ_Advanced.mat'));
                        
        D_spl = cat(3,...
                    data_comb.d_3dX,...
                    data_comb.d_3dY,...
                    data_comb.d_3dZ);
        D_spl = permute(D_spl,[1 3 2]);

        C_abs = cat(2,...
                    data_comb.X_north,...
                    data_comb.Y_east,...
                    data_comb.Z_height);

        L_n = L;

        T_abs = data_comb.time_abs_interf{1};

        save(path_mat,...
             'D_adj',... % Advanced LSQ adjusted Displacement Vectors
             'D_spl',... % Simple LSQ adjusted Displacement Vectors
             'L_adj',... % Advanced LSQ adjusted Line-of-Sight Observations
             'L_n',... % Noisy Simulated Line-of-Sight Observations
             'L_nF',... % Noise-Free Simulated Line-of-Sight Observations (Ground-Truth)
             'C_abs',... % Coordinate of each Point
             'T_abs',... % Timestamps
             'Desc'); % Description of advanced LSQ Adjustement
        
        %% Visualisation
        if create_vis
            N_pt = size(E,1);
            N_coord = size(E,3);
            N_instr = size(L,2);
            N_time = size(L,3);
            vector_scale = 1000;
    
            for fi = 4%1:length(fileName)
                pt_list = 1:5:N_pt;
    
                writerObj = VideoWriter(fullfile(path_out,sprintf('DefVideo_LOS_%s',fileName{fi})));
                writerObj.Quality = 100;
                writerObj.FrameRate = 5;
                open(writerObj);
    
                figure('units','centimeters','position',[0,0,40,20]);
                instr_no = 3;
                for pt_i = pt_list
                    L_adj_plt = squeeze(L_adj{fi}(pt_i,instr_no,:));
                    L_plt = squeeze(L(pt_i,instr_no,:));
                    hold off
                    plot(L_adj_plt,'r','LineWidth',3,'DisplayName','Adjusted');
                    hold on
                    plot(L_plt,'b','LineWidth',3,'DisplayName','Original');
    
                    xlabel('Time');
                    ylabel('Displacement');
                    ylim([min([L_adj{fi}(:);L(:)]),max([L_adj{fi}(:);L(:)])]);
                    legend('Location','southwest')
                    title(sprintf('%s: Observation for Point No. %d',Desc{fi},pt_i));
                    writeVideo(writerObj,getframe(gcf));
                    set(gcf, 'color', [1 1 1]);
                    pause(0.01);
                end
    
                close(gcf)
                close(writerObj);
    
                writerObj = VideoWriter(fullfile(path_out,sprintf('DefVideo_D_%s',fileName{fi})));
                writerObj.Quality = 100;
                writerObj.FrameRate = 5;
                open(writerObj);
    
                figure('units','centimeters','position',[0,0,40,20]);
                title_str = {Desc{fi};'Simple LSQ'};%{'Advanced';'Simple'};
                for s_i = 1:2
                    subplot(1,2,s_i);
                    scatter3(C(:,1),C(:,2),C(:,3),6,'ok','filled');
                    hold on
                end
    
                skip_time = 5; %% Skip Frames for faster Processing <---------------------------------------------
                frame_list = [1,[skip_time:skip_time:N_time]];
                az = linspace(-20,-140,length(frame_list));
                el = linspace(10,5,length(frame_list));
                dx = [min(min(cat(3,data_comb.d_3dX,squeeze(D_adj{fi}(:,1,:)))*vector_scale,[],3),[],2),...
                      max(max(cat(3,data_comb.d_3dX,squeeze(D_adj{fi}(:,1,:)))*vector_scale,[],3),[],2)];
                dy = [min(min(cat(3,data_comb.d_3dY,squeeze(D_adj{fi}(:,2,:)))*vector_scale,[],3),[],2),...
                      max(max(cat(3,data_comb.d_3dY,squeeze(D_adj{fi}(:,2,:)))*vector_scale,[],3),[],2)];
                dz = [min(min(cat(3,data_comb.d_3dZ,squeeze(D_adj{fi}(:,3,:)))*vector_scale,[],3),[],2),...
                      max(max(cat(3,data_comb.d_3dZ,squeeze(D_adj{fi}(:,3,:)))*vector_scale,[],3),[],2)];
                x_limits = [min([C(:,1),C(:,1)]+dx,[],'all'),...
                            max([C(:,1),C(:,1)]+dx,[],'all')];
                y_limits = [min([C(:,2),C(:,2)]+dy,[],'all'),...
                            max([C(:,2),C(:,2)]+dy,[],'all')];
                z_limits = [min([C(:,3),C(:,3)]+dz,[],'all'),...
                            max([C(:,3),C(:,3)]+dz,[],'all')];
    
                for t_i = frame_list
                    C0 = C;
                    C1 = C + D_adj{fi}(:,:,t_i) * vector_scale;
                    C2 = C + [data_comb.d_3dX(:,t_i),...
                              data_comb.d_3dY(:,t_i),...
                              data_comb.d_3dZ(:,t_i)] * vector_scale;
    
                    for pt_i = 1:N_pt
                        Cplt = [C0(pt_i,:);C1(pt_i,:)];
                        Cplt2 = [C0(pt_i,:);C2(pt_i,:)];
                        subplot(1,2,1);
                        plt1(pt_i) = plot3(Cplt(:,1),Cplt(:,2),Cplt(:,3),'-b','LineWidth',2,'DisplayName',title_str{1});
                        subplot(1,2,2);        
                        plt2(pt_i) = plot3(Cplt2(:,1),Cplt2(:,2),Cplt2(:,3),'-m','LineWidth',2,'DisplayName',title_str{2});
                    end
    
                    if t_i == 1
                        ddx = 0.0001*vector_scale;
                        ddyz = 1;
                        for s_i = 1:2
                            ax(s_i) = subplot(1,2,s_i);
                            title(title_str{s_i})
                            view(az(frame_list==t_i),el(frame_list==t_i));
                            axis equal
                            xlim(x_limits + [-ddx ddx]);
                            ylim(y_limits + [-ddx ddx]);
                            zlim(z_limits + [-ddx ddx]);
    
                            grid on
                            box on
                        end
                        set(gcf, 'color', [1 1 1]);
                        Link = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
                        setappdata(gcf, 'StoreTheLink', Link);
                    end
                    for s_i = 1:2
                        view(az(frame_list==t_i),el(frame_list==t_i));
                    end
                    sgtitle(sprintf('%d/%d (Vector Scaling: %dx)',t_i,N_time,vector_scale));
    
                    writeVideo(writerObj,getframe(gcf));
    
                    pause(0.01);
                    if t_i ~= frame_list(end)
                        delete(plt1);
                        delete(plt2);
                    end
    
                end
    
                close(gcf)
                close(writerObj);
            end
        end
        
        fprintf('Finished processing after % 6.1fs\n',toc(t_d0));
    end
end
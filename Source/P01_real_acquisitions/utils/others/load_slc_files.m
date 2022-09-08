function [x_axis, y_axis, complex_data_static, coherence, properties] = load_slc_files(path_dir)
% Script for loading .mat files.

path_mat = dir(fullfile(path_dir,"*.mat"));

mat_j = 0;
set_i = 1;
for mat_i = 1:length(path_mat)
    if ~strcmp([path_mat(mat_i).name],'aoi.mat')
        %% Load File
        fprintf("Loading File %s.\n",path_mat(mat_i).name);
        if ~isfile(fullfile(path_mat(mat_i).folder,path_mat(mat_i).name))
            error('File does not exist. (%s)',...
                  fullfile(path_mat(mat_i).folder,path_mat(mat_i).name));
        end
        mats = load(fullfile(path_mat(mat_i).folder,path_mat(mat_i).name));
        
        %% Loaded File No.
        mat_j = mat_j + 1;
        
        %% Properties
        [~,name,~] = fileparts(path_mat(mat_i).name);
        properties.no_of_file_i(mat_j) = str2double(name(end-4:end));
        properties.length_of_file_i(mat_j)  = size(mats.complex_data_static,2);
    
        if mat_j==1
            file_j(mat_j) = properties.length_of_file_i(mat_j) + 1;
            properties.set_of_file_i(mat_j) = set_i;
        else
            file_j(mat_j) = properties.length_of_file_i(mat_j);
            properties.set_of_file_i(mat_j) = set_i;
            dfile = file_j(mat_j-1)-file_j(mat_j);
            if dfile>1
                set_i = set_i + 1;
            end
        end
            
        [~,n_obs] = size(mats.complex_data_static);
    
        field_names = fields(mats.coh);
        if mat_j==1
            for field_i = 1:length(field_names)
                coherence.(field_names{field_i}) = mats.coh.(field_names{field_i});
            end
        else
            for field_i = 1:length(field_names)
                if ~strcmp(field_names{field_i},'class_id') && ...
                    ~strcmp(field_names{field_i},'class_descr') && ...
                    ~strcmp(field_names{field_i},'coh_1N')
    
                    dims_n = ndims(coherence.(field_names{field_i}));
    
                    if dims_n == 2
                        coherence.(field_names{field_i})(:,end+1) = ...
                            mats.coh.(field_names{field_i});
                    elseif dims_n == 3
                        coherence.(field_names{field_i})(:,:,end+1) = ...
                            mats.coh.(field_names{field_i});
                    end
                elseif strcmp(field_names{field_i},'coh_1N')
                    dims_n = ndims(coherence.(field_names{field_i}));
                    if dims_n == 2
                        [~,nn] = size(mats.coh.(field_names{field_i}));
                        coherence.(field_names{field_i})(:,end+1:end+nn) = ...
                                mats.coh.(field_names{field_i});
                    elseif dims_n == 3
                        [~,~,nn] = size(mats.coh.(field_names{field_i}));
                        coherence.(field_names{field_i})(:,:,end+1:end+nn) = ...
                            mats.coh.(field_names{field_i});
                    end
                
                end
            end
        end
        
        if mat_j==1
            x_axis = mats.x_axis; % Along Track
            y_axis = mats.y_axis; % Across Track
            complex_data_static = mats.complex_data_static; % Static
        else
            complex_data_static(:,end+1:end+n_obs) = mats.complex_data_static; % Static
        end
    end
end

properties.num_of_acquisitions = size(complex_data_static,2);

end


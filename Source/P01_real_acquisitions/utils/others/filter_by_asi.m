function [complex_data,x_axis,y_axis,t_asi,idx_a,idx_b] = filter_by_asi(slc,x_axis,y_axis,id_i,t_asi)

    dim=ndims(slc);
    
    %% Filtering by ASI
    if length(t_asi)==1
        amp = abs(slc);
        amp_sigma = std(amp,1,dim);
        amp_mu = mean(amp,dim);
        asi = 1 - (amp_sigma ./ amp_mu); % Amplitude Stability Index
        filt_asi = asi>=t_asi;
        t_asi = filt_asi; % Overwrite for following loadings.
    else
        filt_asi = t_asi;
    end

    if dim==3
        [idx_a,idx_b] = find(filt_asi); % 2x 1D-Vector
    end
    idx = find(filt_asi); % 2D-Matrix

    if dim==2
        n_idx = length(idx); % Number of Stable Reflectors
    elseif dim==3
        n_idx = length(idx_a); % Number of Stable Reflectors
    end
    n_obs = size(slc,dim); % Number of Observations

    x_axis = x_axis(idx); % Along Track
    y_axis = y_axis(idx); % Across Track

    complex_data = complex(zeros(n_idx,n_obs),0); % Static

    for idx_i = 1:n_idx
        if dim==2
            complex_data(idx_i,:) = slc(idx(idx_i), :);
        elseif dim==3
            complex_data(idx_i,:) = slc(idx_a(idx_i),...
                                        idx_b(idx_i),...
                                        :);
        end
    end
end


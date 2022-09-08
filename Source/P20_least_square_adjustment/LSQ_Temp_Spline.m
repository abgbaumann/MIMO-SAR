function [D_adj, L_adj] = LSQ_Temp_Spline(L,E,S,cov,T)
% This functions takes as an input two matrices containing Line-of-Sight
% displacements and unit vectors.
%
% -------------------------------------------------------------------------
%
% L:            [p x i x t] - LOS-Displacements over time for each instrument
% E:            [p x i x c] - Unit vector (3D) for each point to each
%                             instrument
% S:            [p x i]     - Standard deviation of each point for each
%                             instrument
% cov:          [c x 2]     - Processing values for parametric function
%                             (e.g. [1,20;1,20;1,20])$
% T:            [t x 1]     - Vector of timestamps [optional]
% 
% - with p being the number of points in the point cloud
% - with i being the number of instruments observing the points
% - with t being the number of acquisitions (i.e. time)
% - with c being the number of coordinate components (i.e. 3 for X,Y,Z)
%
% -------------------------------------------------------------------------
%
% D_adj:        [p x c x t] - 3D-Displacement-Vector for each point and
%                             time.
% L_adj:        [p x i x t] - LOS-Displacements after adjustment
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang & Jemil Avers Butt, ETH Zürich, GSEG (26th July 2022)
%

plt_fig = 0; % Plot Figures
plt_ani = 0; % Plot Animation - plt_fig needs to be activated too; script becomes slow if activated.

step_i = 0;
tic0 = tic;

N_pt = size(E,1);
N_coord = size(E,3);
N_instr = size(L,2);
N_time = size(L,3);

if ~exist('T','var')
    T = 1:N_time;
end

%% Observations (LOS)
step_i = step_i+1;
fprintf("Step % 2d (% 6.1fs): Observation-Vector\n",step_i,round(toc(tic0),1));        
for pt_i = 1:N_pt
    b_tmp = squeeze(L(pt_i,:,:));
    b(:,pt_i) = b_tmp(:); 
end

%% Matrices
step_i = step_i+1;
fprintf("Step % 2d (% 6.1fs): Functional Model\n",step_i,round(toc(tic0),1)); 
% Cov funs
s_x = cov(1,1);    d_x = cov(1,2);
s_y = cov(2,1);    d_y = cov(2,2);
s_z = cov(3,1);    d_z = cov(3,2);

k_fun_x=@(s,t) s_x^2*exp(-(s-t).^2/(d_x).^2);
k_fun_y=@(s,t) s_y^2*exp(-(s-t).^2/(d_y).^2);
k_fun_z=@(s,t) s_z^2*exp(-(s-t).^2/(d_z).^2);

% untransformed cov_mats
K_x=zeros(N_time,N_time);
K_y=zeros(N_time,N_time);
K_z=zeros(N_time,N_time);

for k=1:N_time
    for l=1:N_time
        K_x(k,l)=k_fun_x(T(k),T(l));
        K_y(k,l)=k_fun_y(T(k),T(l));
        K_z(k,l)=k_fun_z(T(k),T(l));
    end
end
 
% transformed_cov_mat
K_ij=blkdiag(K_x,K_y,K_z);

%% Matrix of linear functionals
step_i = step_i+1;
fprintf("Step % 2d (% 6.1fs): Design-Matrix\n",step_i,round(toc(tic0),1));  

A = zeros(N_time*N_instr,N_time*N_coord,N_pt);
for pt_i = 1:N_pt
    A_x_tmp=squeeze(E(pt_i,:,1))';
    A_y_tmp=squeeze(E(pt_i,:,2))';
    A_z_tmp=squeeze(E(pt_i,:,3))';
    
    A_slice=[A_x_tmp, zeros(N_instr,N_time-1), A_y_tmp,zeros(N_instr,N_time-1), A_z_tmp,zeros(N_instr,N_time-1)];    
    for k=1:N_time
        A(1+(k-1)*N_instr:1+k*N_instr-1,:,pt_i)=circshift(A_slice, [0,k-1]);
    end
end

%% Optimal estimation
step_i = step_i+1;
fprintf("Step % 2d (% 6.1fs): Optimal Estimation\n",step_i,round(toc(tic0),1));  

D_adj = zeros(N_pt,N_coord*N_time);
L_adj = zeros(N_pt,N_instr*N_time);

fbar = waitbar(0,'Optimal Estimation...');
for pt_i = 1:N_pt
    % Setup matrices

    % Define noise variance = if higher, then smoother (degree of
    % noninterpolation)
    s = repmat(S(pt_i,:),1,N_time);
    P = diag(s.^2);

    K_ij_transformed = A(:,:,pt_i)*K_ij*A(:,:,pt_i)';
    K_t_transformed = A(:,:,pt_i)*K_ij;

    % Solve estimation equation
    %t0=tic;lambda_vec1= pinv(K_ij_transformed+P)*b(:,pt_i);toc(t0)
    %t0=tic;lambda_vec2= lsqminnorm(K_ij_transformed+P,b(:,pt_i));toc(t0)
    %t0=tic;lambda_vec3= (K_ij_transformed+P)\b(:,pt_i);toc(t0)
    %t0=tic;lambda_vec4= linsolve(K_ij_transformed+P,b(:,pt_i));toc(t0)
    lambda_vec = pinv(K_ij_transformed+P)*b(:,pt_i);
    
    f_hat = lambda_vec'*K_t_transformed;

    D_adj(pt_i,:) = f_hat;
    
    L_adj(pt_i,:) = A(:,:,pt_i)*f_hat';
    
    waitbar(pt_i/N_pt,fbar,'Optimal Estimation...');    
end
close(fbar);

%% Reshape for return
D_adj = reshape(D_adj,N_pt,N_time,N_coord);
D_adj = permute(D_adj,[1 3 2]);
L_adj = reshape(L_adj,N_pt,N_instr,N_time);

%% Plotting 1
if plt_fig
    if plt_ani
        pt_list = 1:N_pt;
    else
        pt_list = floor(N_pt/2);
    end

    
    for n_i = 1:3
        figure('units','centimeters','position',[0,0,40,20]);
    for pt_i = pt_list
        L_adj_plt = squeeze(L_adj(pt_i,n_i,:));
        L_plt = squeeze(L(pt_i,n_i,:));
        hold off
        plot(T,L_adj_plt,'r','DisplayName','Adjusted');
        hold on
        plot(T,L_plt,'b','DisplayName','Original');

        xlabel('Time');
        ylabel('Displacement');
        ylim([min([L_adj(:);L(:)]),max([L_adj(:);L(:)])]);
        legend('Location','best')
        title(sprintf('Observation for Point No. %d',pt_i));
        pause(0.01);
    end
    end
end

end

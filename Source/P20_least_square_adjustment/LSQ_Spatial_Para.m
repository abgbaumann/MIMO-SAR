function [D_adj, L_adj] = LSQ_Spatial_Para(L,E,S,fmod)
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
% fmod:         cell        - Defines parametric function in the first cell
%                             (e.g. 'poly02': i.e. y = a + b*x + c*x^2) 
%                             and the projection vector(s) in the second 
%                             cell.
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

%% Functional Model
if strfind(fmod{1},'poly')
    fdeg = str2double(fmod{1}(5:end));
else
    error('Unknown Function');
end

if size(fmod{2},1)>1
    C = fmod{2};
else
    C = repmat(fmod{2},size(E,1),1);
end

%% Observations (LOS)
step_i = step_i+1;
fprintf("Step % 2d (% 6.1fs): Observation-Vector\n",step_i,round(toc(tic0),1));        
for t_i = 1:N_time
    b_tmp = squeeze(L(:,:,t_i))';
    b(:,t_i) = b_tmp(:); 
end

%% Design-Matrix for Projection with Parametric Conditions from LOS to 3D Displacement
N_para = 3*(fdeg+1);
h=zeros(N_coord,N_para*N_coord);
for p_i = 1:N_para
    for c_i = 1:N_coord
        i = (c_i-1)*N_para + p_i; 
        h(c_i,i) = 1;
    end
end

fdeg_vec = repmat([0:fdeg],N_coord,1);
fdeg_vec = fdeg_vec(:);
H = repmat(h,N_pt,1);
for para_i = 1:size(H,2)
    d_i = rem(para_i,N_para);
    if rem(d_i,N_coord) == 0
        d_i = N_coord;
    else
        d_i = rem(d_i,N_coord);
    end
    if rem(para_i,N_para) == 0
        f_i = N_para;
    else
        f_i = rem(para_i,N_para);
    end
    H(H(:,para_i)==1,para_i) = C(:,d_i).^fdeg_vec(f_i);
end

%% Design-Matrix for Backprojection from 3D Dispalcement to LOS
for pt_i = 1:N_pt
    R_tmp{pt_i} = squeeze(E(pt_i,:,:));
end
R = blkdiag(R_tmp{:});

clearvars('R_tmp','pt_i','d_i','hi','t_i','b_tmp');

%% Processing
A = R*H;

if ~(size(S,1)==N_instr)
    S=S';
end

P = diag(1./S(:).^2);

x = zeros(N_para*N_coord,N_time);

L_adj = zeros(N_instr*N_pt,N_time);
D_adj = zeros(N_coord*N_pt,N_time);

APAinv = pinv(A.'*P*A);
for t_i = 1:N_time
    l = b(:,t_i);
    x(:,t_i) = (APAinv*A.'*P*l)';

    L_adj(:,t_i) = A*x(:,t_i);
    X = H*x(:,t_i);

    D_adj(:,t_i) = X;
end

%% Reshape for return
D_adj = reshape(D_adj,N_coord,N_pt,N_time);
D_adj = permute(D_adj,[2,1,3]);
L_adj = reshape(L_adj,N_instr,N_pt,N_time);
L_adj = permute(L_adj,[2 1 3]);

%% Plotting 1
if plt_fig
    if plt_ani
        pt_list = 1:N_pt;
    else
        pt_list = floor(N_pt/2);
    end

    figure('units','centimeters','position',[0,0,40,20]);
    for pt_i = pt_list
        L_adj_plt = squeeze(L_adj(pt_i,2,:));
        L_plt = squeeze(L(pt_i,2,:));
        hold off
        plot(L_adj_plt,'r','DisplayName','Adjusted');
        hold on
        plot(L_plt,'b','DisplayName','Original');

        xlabel('Time');
        ylabel('Displacement');
        ylim([min([L_adj(:);L(:)]),max([L_adj(:);L(:)])]);
        legend('Location','best')
        title(sprintf('Observation for Point No. %d',pt_i));
        pause(0.01);
    end
end

end


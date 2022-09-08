function [D_adj, L_adj] = LSQ_Temp_Para(L,E,S,fmod,T)
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
% fmod:         cell        - Defines parametric function (e.g. 'poly02': 
%                             y = a + b*x + c*x^2)
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
% by Andreas Baumann-Ouyang & Jemil Avers Butt, ETH Zürich, GSEG (05 May 2022)
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

%% Design-Matrix for Projection with Parametric Conditions from LOS to 3D Displacement
% Functional model
step_i = step_i+1;
fprintf("Step % 2d (% 6.1fs): Functional Model\n",step_i,round(toc(tic0),1));  
if strfind(fmod{1},'poly')
    fdeg = str2double(fmod{1}(5:end));
else
    error('Unknown Function');
end

N_para = (fdeg+1);
g=zeros(N_coord,N_para*N_coord);
for p_i = 1:N_para
    for c_i = 1:N_coord
        i = (c_i-1)*N_para + p_i; 
        g(c_i,i) = 1;
    end
end

fdeg_vec = repmat([0:fdeg]',N_coord,1);
fdeg_vec = fdeg_vec(:);
G = repmat(g,N_time,1);
for para_i = 1:size(G,2)
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

    G(G(:,para_i)==1,para_i) = T.^fdeg_vec(f_i);
end

%% Design-Matrix for Backprojection from 3D Dispalcement to LOS
step_i = step_i+1;
fprintf("Step % 2d (% 6.1fs): Design-Matrix\n",step_i,round(toc(tic0),1));  
Q = zeros(N_time*N_instr,N_time*N_coord,N_pt);
for pt_i = 1:N_pt
    Q_tmp = squeeze(E(pt_i,:,:));
    
    Q_tmp2 = repmat(Q_tmp, 1, N_time);
    Q_tmp3 = mat2cell(Q_tmp2, N_instr, repmat(N_coord,1,N_time));
    Q(:,:,pt_i) = blkdiag(Q_tmp3{:});
end
clearvars('Q_tmp','pt_i','d_i','gi', 'b_tmp');

%% Processing
step_i = step_i+1;
fprintf("Step % 2d (% 6.1fs): Optimal Estimation\n",step_i,round(toc(tic0),1));  
D_adj = zeros(N_pt,N_coord*N_time);
L_adj = zeros(N_pt,N_instr*N_time);

for pt_i = 1:N_pt
    A = Q(:,:,pt_i)*G;
    s = repmat(S(pt_i,:),1,N_time);
    P = diag(1./s.^2); %P = diag(s0^2./s.^2);
    
    l = b(:,pt_i);
    
    APAinv = pinv(A.'*P*A);
    x = (APAinv*A.'*P*l)';

    L_adj(pt_i,:) = A*x';
    X = G*x';

    D_adj(pt_i,:) = X;
end

%% Reshape for return
D_adj = reshape(D_adj,N_pt,N_coord,N_time);
L_adj = reshape(L_adj,N_pt,N_instr,N_time);

%% Plotting 1
if plt_fig
    if plt_ani
        pt_list = 1:N_pt;
    else
        pt_list = floor(N_pt/2);
    end

    figure('units','centimeters','position',[0,0,40,20]);
    instr_no = 1;
    for pt_i = pt_list
        L_adj_plt = squeeze(L_adj(pt_i,instr_no,:));
        L_plt = squeeze(L(pt_i,instr_no,:));
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


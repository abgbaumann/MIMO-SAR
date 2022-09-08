function [X_N,Y_E,Z_H,dx,dy,dz] = random_vector_field_deformations(E_0,E_1,E_N,N_0,N_1,N_N,H_0,H_1,H_N)
% This functions simulates a random vector field deformation.
%
% -------------------------------------------------------------------------
% by Jemil Avers Butt, ETH Zürich, GSEG (8th September 2022)
%

%
% E_0 = 0;
% E_1 = 1;
% E_N = 20;
% N_0 = 0;
% N_1 = 1;
% N_N = 20;
% H_0 = 0;
% H_1 = 1;
% H_N = 100;
% [X_N,Y_E,Z_H,dx,dy,dz] = random_vector_field_deformations(E_0,E_1,E_N,N_0,N_1,N_N,H_0,H_1,H_N)
%

plt_fig = 0;
plt_fig3d = 0;

% Define some  basic parameters
approximate='yes' ;   % yes or no for truncation after explained_var of variance
explained_var=0.99;

sigma_0_h=1;

e_vector=linspace(E_0,E_1,E_N);
n_vector=linspace(N_0,N_1,N_N);
h_vector=linspace(H_0,H_1,H_N);

sigma_e = @(x1,x2) exp(-abs(((x1-x2)/0.3).^2));  % squared exponential
sigma_n = @(y1,y2) exp(-abs(((y1-y2)/0.3).^2));  % squared exponential
sigma_h = @(z1,z2) sigma_0_h^2*exp(-abs(((z1-z2)/0.15).^2));  % squared exponential

R_e=bsxfun(sigma_e,e_vector,e_vector');
[U_e, S_e, ~]=svd(R_e);
e_e_mat=U_e;   % one column = one eigenfunction evaluated along x

R_n=bsxfun(sigma_n,n_vector,n_vector');
[U_n, S_n, ~]=svd(R_n);
e_n_mat=U_n;   % one column = one eigenfunction evaluated along y

R_h=bsxfun(sigma_h,h_vector,h_vector');
[U_h, S_h, ~]=svd(R_h);
e_h_mat=U_h;   % one column = one eigenfunction evaluated along y

for ww=1:3

    % Now start simulating:
    % 1. Make N(0,1) RV's
    q=normrnd(0,1,[E_N*N_N*H_N,1]);

    % 2. Calculate and sort the eigenvalues of the tensor product kernel
    LAMBDA=zeros(E_N*N_N*H_N,1);  %  eigenvalues,
    Ind=zeros(E_N*N_N*H_N,3);     %  i, j index

    for i=1:E_N
        for j=1:N_N
           for k=1:H_N
            LAMBDA((i-1)*E_N*N_N+(j-1)*N_N+k,1)=S_e(i,i)*S_n(j,j)*S_h(k,k);
            Ind((i-1)*E_N*N_N+(j-1)*N_N+k,:)=[i,j,k];
           end
        end
    end

    [LAMBDA_sort , index]=sortrows(LAMBDA,[-1]);% Sort descending eigenvalues
    Ind_sort=Ind(index,:);                    % Remember sorting for eigenfunctions 

    tot_var=sum(LAMBDA_sort);
    acc_lambda=cumsum(LAMBDA_sort);
    trunc_ind=find(acc_lambda/tot_var>=explained_var,1);

    Acc=zeros(E_N,N_N,H_N);
    delta_st=zeros(E_N,N_N,H_N);

    if strcmp(approximate,'yes')==1
    for k=1:trunc_ind
        delta=q(k)*sqrt(LAMBDA_sort(k))*e_n_mat(:,Ind_sort(k,1))*e_e_mat(:,Ind_sort(k,2))';  % x,y change because of illustration
        for i=1:H_N
            delta_st(:,:,i)=e_h_mat(i,Ind_sort(k,3))*delta;
        end  
        Acc=Acc+ delta_st;
    end
    else
        for k=1:E_N*N_N*H_N
        delta=q(k)*sqrt(LAMBDA_sort(k))*e_n_mat(:,Ind_sort(k,1))*e_e_mat(:,Ind_sort(k,2))';  % x,y change because of illustration
        for i=1:H_N
            delta_st(:,:,i)=e_h_mat(i,Ind_sort(k,3))*delta;
        end  
        Acc=Acc+ delta_st;
        end
    end

    if ww==1
        e_vals_tor=Acc;
    elseif ww==2
        n_vals_tor=Acc;
    elseif ww==3
        h_vals_tor=Acc;
    end

end

[ee,nn]=meshgrid(e_vector',n_vector');
hh=zeros(size(ee,1),size(ee,2));

if plt_fig
    xlimits_val = [min(ee(:))-0.1,max(ee(:))+0.1];
    ylimits_val = [min(nn(:))-0.1,max(nn(:))+0.1];
    zlimits_val = [min(hh(:))-0.1,max(hh(:))+0.1];
    
    opengl software

    % Initialize the Videowriter
    if plt_fig3d
        str_add = '3d';
    else
        str_add = '2d';
    end
    writerObj = VideoWriter(sprintf('Random_deformation_field_%s',str_add));
    writerObj.Quality = 100;
    writerObj.FrameRate = 30;
    open(writerObj);

    hFig=figure('units','normalized','outerposition',[0 0 1 1]);
    for k=2:H_N
        e_vals_mat=squeeze(e_vals_tor(:,:,k));
        n_vals_mat=squeeze(n_vals_tor(:,:,k));
        h_vals_mat=squeeze(h_vals_tor(:,:,k));
        
        e_vals_mat_pre=squeeze(e_vals_tor(:,:,k-1));
        n_vals_mat_pre=squeeze(n_vals_tor(:,:,k-1));
        h_vals_mat_pre=squeeze(h_vals_tor(:,:,k-1));
        
        de_vals_mat = e_vals_mat_pre - e_vals_mat;
        dn_vals_mat = n_vals_mat_pre - n_vals_mat;
        dh_vals_mat = h_vals_mat_pre - h_vals_mat;
        
        subplot(1,2,1);
        if plt_fig3d
            quiver3(ee(:),nn(:),hh(:),e_vals_mat(:),n_vals_mat(:),h_vals_mat(:));
            zlim(zlimits_val);
        else
            quiver(ee(:),nn(:),e_vals_mat(:),n_vals_mat(:));
        end
        xlim(xlimits_val);
        ylim(ylimits_val);
        xlabel('E');ylabel('N');zlabel('H');
        title(sprintf('Absolute Displacement [t_{%03d}]',k));
        
        subplot(1,2,2);
        
        if plt_fig3d
            quiver3(ee(:),nn(:),hh(:),de_vals_mat(:),dn_vals_mat(:),dh_vals_mat(:));
            zlim(zlimits_val);
        else
            quiver(ee(:),nn(:),de_vals_mat(:),dn_vals_mat(:));
        end
        xlim(xlimits_val);
        ylim(ylimits_val);
        title(sprintf('Relative Displacement [t_{%03d} - t_{%03d}]',k-1,k));
        
        pause(0.01)
        set(gcf, 'color', [1 1 1]);
        
        writeVideo(writerObj,getframe(gcf));
        
    end

    close(gcf)
    close(writerObj);
end

Y_E = ee(:)';
X_N = nn(:)';
Z_H = hh(:)';

dy = e_vals_tor;
dy = reshape(dy,size(Y_E,2),[])';

dx = n_vals_tor;
dx = reshape(dx,size(X_N,2),[])';

dz = h_vals_tor;
dz = reshape(dz,size(Z_H,2),[])';

end



















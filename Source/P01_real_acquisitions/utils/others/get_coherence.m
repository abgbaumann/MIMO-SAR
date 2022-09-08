function[Coherence_mat]=get_coherence(SLC_1,SLC_2,box_size)
%
% This function calculates the Coherence of two SLC's. for averaging
% purposes the box_size defines the moving window used for estimation.
% SLC_1 and SLC_2 need to have the same dimension.
%
% For this do the following:
%
% 1. Definitions and imports
% 2. Moving window filter  
% 
% The inputs are:
% Name                         meaning                           format
% SLC_1             SLC image (complex matrix)                 [n_row,n_col]
% SLC_2             SLC image (complex matrix)                 [n_row,n_col]
% box_size          Window size for averaging                  [2,1]
%
% The outputs are:
% Name                         meaning                           format
% Coherence mat     Matrix of real coherence values            [n_row,n_col]
 
% 1. Definitions and imports ---------------------------------------------
 
SLC_1_1=SLC_1.*conj(SLC_1);
SLC_2_2=SLC_2.*conj(SLC_2);
SLC_1_2=SLC_1.*conj(SLC_2);
 
 
% 2. Moving window filter  
 
SLC_1_1_avg=imfilter(SLC_1_1, ones(box_size), 'same');
SLC_2_2_avg=imfilter(SLC_2_2, ones(box_size), 'same');
SLC_1_2_avg=imfilter(SLC_1_2, ones(box_size), 'same');
 
Coherence_mat=abs(SLC_1_2_avg)./(sqrt(abs(SLC_1_1_avg)).*(sqrt(abs(SLC_2_2_avg))));
 
 
end

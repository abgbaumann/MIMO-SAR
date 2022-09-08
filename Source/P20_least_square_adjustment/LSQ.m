function [x] = LSQ(L,A,P)
% This functions takes as an input the observation L, coefficient matrix A,
% and weight matrix P and performs a least square adjustment.
%
% -------------------------------------------------------------------------
%
% Observation Equation:         v = A * x - L
%
% Normal Equation:              N = A^T * P * A
%
% h-Vector:                     h = A^T * P * L
%
% Vector of Unknown Parameters: x = N^(-1) * h
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich, GSEG (07th March 2022)
%

x = ((A.'*P*A)^(-1)*A.'*P*L)'; % Solution of the Standard Equation

end


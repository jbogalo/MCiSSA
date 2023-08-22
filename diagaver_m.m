function y = diagaver_m(Y,M)
% DIAGAVER_M - Diagonal averaging for Multivariate Singular Spectrum Analysis.
%
% This MATLAB function transforms the matrix, Y, into the multivariate
% time series, y, with M time series by columns, by diagonal averaging.
% This entails averaging the vector elements of Y over its antidiagonals.
%
% Syntax:     y = diagaver_m(Y,M)
%
% Input arguments:
% Y:    Matrix associated with a group of frequencies.
% M:    Number of time series.
%
% Output arguments:
% y:    Matrix with M columns containing the reconstructed components.

% -------------------------------------------------------
% Realignment
% -------------------------------------------------------
[LL, NN] = size(Y);
if  mod(LL,M)||LL/M>NN
    Y = Y';
end

% -------------------------------------------------------
% Dimensions
% -------------------------------------------------------
[LL, N] = size(Y);
L = LL/M;
T = N+L-1;

% -------------------------------------------------------
% Diagonal averaging
% -------------------------------------------------------
y = zeros(T,M);
for i=1:M
    v = zeros(1,M);
    v(i) = 1;
    P = kron(eye(L),v);
    y(:,i) = diagaver(P*Y);
end

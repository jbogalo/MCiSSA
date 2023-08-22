function [V, D, Zs, Z] = mcissa(x,L,varargin)
% M-CiSSA - Multivariate Circulant Singular Spectrum Analysis.
%
% This MATLAB function, using Multivariate Circulant Singular Spectrum
% Analysis, returns the elementary reconstructed subcomponents and
% components by frequency, Z and Z respectively. Also, it computes the
% eigenvalues, D, and real eigenvectors, V, of block circulant matrix
% obtained with the second moments of the input vector time series, x, 
% given a window length, L.
%
% Syntax:     [V, D, Zs, Z] = mssa_jb(x,L)
%             [V, D, Zs, Z] = mssa_jb(x,L,H)
%
% Input arguments:
% x: Matrix containing the M time series original data by columns.
% L: Window length.
% H: Optional. Related to the characteristics of the time series.
%    H=0 Autoregressive extension (default). It is indicated for stationary
%        and stochastic trend time series as well.
%    H=1 Mirroring. It can be used with stationary time series and works
%        well for AM-FM signals.
%    H=2 No extension. It is suitable for deterministic time series.
%
% Output arguments:
% V:  Real eigenvectors of block circulant matrix.
% D:  Cell array with the L eigenvalues matrices of the cross 
%     spectral density matrices.
% Zs: Cell array with the M matrices containing the elementary
%     reconstructed subcomponents by frequency.
% Z:  Cell array with the M matrices containing the elementary
%     reconstructed components by frequency.

% -------------------------------------------------------
% Checking the input arguments
% -------------------------------------------------------
% Dimensions
[T, M] = size(x);
N = T-L+1;
if L>N
    error(' ***  The window length must be less than T/2  *** ');
end
LM = L*M;

% Type of extension depending on H
if nargin>2
    switch varargin{1}
        case 1
            H = T;
        case 2
            H = 0;
        otherwise
            H = L;
    end
else
    H = L;
end

% Number of symmetryc frequency pairs around 1/2
if mod(L,2)
    nf2 = (L+1)/2-1;
else
    nf2 = L/2-1;
end

% Number of frequencies <= 1/2
nft = nf2+abs(mod(L,2)-2);

% -------------------------------------------------------
% Trajectory matrix
% -------------------------------------------------------
% Extended series
xe = zeros(T+2*H,M);
for i=1:M
    xe(:,i) = extend(x(:,i),H);
end

% Trajectory matrix
X = zeros(LM,N+2*H);
for j=1:L
    X((j-1)*M+1:j*M,:) = xe(j:end-L+j,:)';
end

% -------------------------------------------------------
% Decomposition
% -------------------------------------------------------
% Cross-covariance matrix function
Gam = zeros(M,M,L);
for m=0:L-1
    Gam(:,:,m+1)=(x(1:T-m,:)-kron(mean(x),ones(T-m,1)))'*(x(1+m:T,:)-kron(mean(x),ones(T-m,1)))/(T-m);
end

% Block Toeplitz cross-covariance matrix S and equivalent block circulant matrix C
S = kron(eye(L),Gam(:,:,1)); C = S;
for j=1:L
    for k=j+1:L
        m = abs(j-k);
        S((j-1)*M+1:j*M,(k-1)*M+1:k*M) = Gam(:,:,m+1);
        S((k-1)*M+1:k*M,(j-1)*M+1:j*M) = Gam(:,:,m+1)';
        C((j-1)*M+1:j*M,(k-1)*M+1:k*M) = ((L-m)/L)*Gam(:,:,m+1)+(m/L)*Gam(:,:,L-m+1)';
        C((k-1)*M+1:k*M,(j-1)*M+1:j*M) = ((L-m)/L)*Gam(:,:,m+1)'+(m/L)*Gam(:,:,L-m+1);
    end
end
clear Gam

% Fourier unitary matrix
U = dftmtx(L)/sqrt(L);

% Cross spectral density matrices
F = cell(L,1);
Fc = sqrt(L)*kron(U,eye(M))'*C(:,1:M);
for k=1:L
    F{k} = Fc((k-1)*M+1:k*M,:);
end

% Diagonalization of cross spectral density matrices
E = cell(L,1); D = cell(L,1);
for k=1:L
    [E{k}, D{k}] = eigs(F{k},M);
    D{k} = abs(D{k});
end

% Eigenvectors of block circulant matrix (unitary base)
V = kron(U,eye(M))*blkdiag(E{:});

% Real eigenvectors (orthonormal base)
V(:,1:M) = real(V(:,1:M));
for k=1:nf2
    v_l = V(:,k*M+1:(k+1)*M);
    V(:,k*M+1:(k+1)*M) = sqrt(2)*real(v_l);
    V(:,(L-k)*M+1:(L-k+1)*M) = sqrt(2)*imag(v_l);
end
if ~mod(L,2)
    V(:,M*(nft-1)+1:M*nft) = real(V(:,M*(nft-1)+1:M*nft));
end

% Principal components
W = V'*X;

% -------------------------------------------------------
% Reconstruction
% -------------------------------------------------------
% Elementary reconstructed subcomponents
Rs = cell(M,1);
for j=1:M
    Rs{j} = zeros(T+2*H,M,L);
end
for k=1:L
    for m=1:M
        y = diagaver_m(V(:,(k-1)*M+m)*W((k-1)*M+m,:),M);
        for j=1:M
            Rs{j}(:,m,k) = y(:,j);
        end
    end
end

% Elementary reconstructed components
R = cell(M,1);
for j=1:M
    R{j} = sum(Rs{j},2);
end

% -------------------------------------------------------
% Grouping by frequency
% -------------------------------------------------------
% Elementary reconstructed subcomponents by frequency
Zs = cell(M,1);
for j=1:M
    Zs{j} = zeros(T+2*H,M,nft);
    Zs{j}(:,:,1) = Rs{j}(:,:,1);
    for k=1:nf2
        Zs{j}(:,:,k+1) = Rs{j}(:,:,k+1)+Rs{j}(:,:,L+1-k);
    end
    if ~mod(L,2)
        Zs{j}(:,:,nft) = Rs{j}(:,:,nft);
    end
    Zs{j} = Zs{j}(H+1:end-H,:,:);
end    

% Elementary reconstructed components by frequency
Z = cell(M,1);
for j=1:M
    Z{j} = sum(Zs{j},2);
end

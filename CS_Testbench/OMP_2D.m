function Z_hat = OMP_2D(A,Y)  
% INPUT:
%   A = m x n sampling matrix
%   Y = m x m sample
%   k = (optional) sparsity level
% OUTPUT:
%   Z_hat = n x n reconstruction of the ideal signal Z

% OTHER VARIABLES:
%   R = m x m residual measurements
%   (i, j) = set of the coordinates of atoms that are allowed to be
%   selected in the future, i for row indices and j for column indices

% idx = it, nstar = bold i, u_hat = x_val

% Initialization
maxiter = size(Y,2);
tol = 1e-6;
k = 10;                      % sparsity level
[m, n] = size(A);            % m:dim of signal, n:#atoms in dictionary

R_old = Y;                  % residual of y
Z_hat = zeros(n,n);         % coefficient (output)
An = [];                    % support set
nstar = [];                 % selected support

for t=1:maxiter

    [~,idx] = max(abs(A'*R_old));   % scan for best candidate
    nstar(t,:) = ind2sub(size(Y),idx);
    An = [An A(:,idx)];             % augment desired support set
    x_val = pinv(An)*y;             % best signal estimate so far
    R_new = y - An*x_val;           % update signal residue
    
    Z_hat(nstar) = x_val;
    
    % Stopping criteria
    if ((norm(A*Z_hat - y) < tol) || (norm(R_new) > norm(r_old)))
        break;
    end
    
    r_old = R_new;
end

end
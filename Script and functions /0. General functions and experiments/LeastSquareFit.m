%% THE FUNCTION LeastSquareFit COMPUTE THE FITTING ERROR OF THE LEAST SQUARE 
%% APPROXIMATED B-SPLINE DEFINED ON A GIVEN KNOT VECTOR.
%
% It follows the algorithm 2: Least Square Fit from Kang et al. (2015)
% The function works with non-parametric and parametric curves (open).

function err = LeastSquareFit(U_inner, p, x, y, a, b)
% INPUTS:
% - U_inner: vector of inner knots
% - p: spline degree
% - x: parameters (Nx1)
% - y: data values -> (Nx1) for scalar functions or (NxD) for parametric curves
% - a,b: domain end-points
% ---------------------------------------------
% OUTPUT:
% - err: fitting error
% ---------------------------------------------


% Force x to be a column vector
x = x(:);

% y can be Nx1 (scalar) OR NxD (parametric). Do NOT vectorize y.

% If y is a vector => force it to be a column vec.
if (size(y,1)==1 || size(y,2) == 1)      
    y = y(:);                
end

% Construct the full (clamped) knot vector:
U_full = [a*ones(1, p+1), sort(U_inner(:).'), b*ones(1, p+1)];

n_dim = length(U_full) - p - 1;   % number of basis functions (= #control points)
N = length(x);  % Number of data points

% Build collocation matrix A (N x n_dim):
A = sparse(N, n_dim);

% Loop over the data points:
for k = 1:N
    ik = findspan(p, x(k), n_dim-1, U_full);        % Find the knot span of U_init where x(k) lives
    bv = nonvanishing_basis(ik, p, x(k), U_full);   % Find the non vanishing basis fcts (evaluated in x(k))
    A(k, (ik-p):ik) = bv;        % Save these values of the fcts in the k th row of A and corresponding columns (i_k-p ... i_k)
end

% Least squares solve:
% - scalar: C is (n_dim x 1)
% - parametric: C is (n_dim x D)
C = A \ y;

% Residuals:
R = A*C - y;                   % Nx1 or NxD

% Mean squared residual per sample:
% scalar: sum(R.^2)/N
% parametric: sum over all components (Frobenius norm)
err = (norm(R,'fro')^2) / N;
end


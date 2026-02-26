%% THE FUNCTION LeastSquareFit COMPUTE THE FITTING ERROR OF THE LEAST SQUARE 
%% APPROXIMATED B-SPLINE DEFINED ON A GIVEN KNOT VECTOR.
%
% It follows the algorithm 2: Least Square Fit from Kang et al. (2015)
% The function works with closed periodic curves.

function err = LeastSquareFit_closed(U_inner, p, x, y, a, b)
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

% Construct the full (periodic) knot vector:
K = [a, U_inner, b]; % Add the end-points
dK = diff(K);        % Compute the length of each sub-interval

% Periodic extention of the vector:
left_ext  = a - cumsum(dK(end:-1:end-p+1));
left_ext  = left_ext(end:-1:1);
right_ext = b + cumsum(dK(1:p));

U_full = [left_ext, K, right_ext];      % Full periodic knot vector

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

M = n_dim - p;  % Number of distinct CP
A_final = A(:, 1:M); % Initialize A_final as the first M columns of A

% Loop over the first p columns of A_final
for r = 1:p
    A_final(:, r) = A_final(:, r) + A(:, M + r);
end

% LS solution for the distinct control points:
% Original problem: min_Q_final ||A_final*Q_final - y||
% But if A is ill conditioned o close to be rank deficient NO
% => We use Tikonov:
% min_Q_final ||A_final*Q_final - y||+ epsilon ||Q_final||
% => we have to solve the linear system: (At* A + esp Id)*Q = At * y
% => Subsitute Q_final = A_final \ y with:
epsilon = 1e-5;         % Regularization parameter
ATA = A_final' * A_final;
ATb = A_final' * y;
Q_final = (ATA + epsilon * speye(size(ATA))) \ ATb;

C = [Q_final; Q_final(1:p,:)];   % full periodic CP

% Residuals:
R = A*C - y;                   % Nx1 or NxD

% Mean squared residual per sample:
% scalar: sum(R.^2)/N
% parametric: sum over all components (Frobenius norm)
err = (norm(R,'fro')^2) / N;
end


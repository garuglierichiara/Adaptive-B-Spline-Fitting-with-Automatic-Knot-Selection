%% THE FUNCTION LeastSquareFit_Surf COMPUTES THE Fitting LS ERROR IN THE BIVARIATE CASE BASED ON KANG 2015

function [err] = LeastSquareFit_Surf(U_inner, V_inner, p, q, u_vec, v_vec, Sx, Sy, Sz, ua, ub, va, vb)
% INPUTS:
% - U_inner: knots in the u direction
% - V_inner: knots in the v direction
% - p: degree in the u direction
% - q: degree in the v direction
% - u_vec: evaluation points in the u direction
% - v_vec: evaluation points in the v direction
% - Sx, Sy, Sz: surface components
% - direction: 'u' if we are optimizing inner_u, 'v' if we are optimizing inner_v
% - ua, ub, va, vb: domain end-points
% --------------------------------------------
% OUTPUT:
% - err: fitting error
% --------------------------------------------

% Knot vectors clamped between 0 e 1
U_full = [ua*ones(1,p+1), sort(U_inner(:).'), ub*ones(1,p+1)];
V_full = [va*ones(1,q+1), sort(V_inner(:).'), vb*ones(1,q+1)];

n_dimU = length(U_full) - p - 1;
n_dimV = length(V_full) - q - 1;
[n_rows, n_cols] = size(Sz); 
N_tot = n_rows * n_cols;

Target = [Sx(:), Sy(:), Sz(:)];
A = sparse(N_tot, n_dimU * n_dimV);

% Costrustruction of the matrix A using the tensor product
for k = 1:n_cols
    u = u_vec(k);
    iu = findspan(p, u, n_dimU-1, U_full);
    Nu = nonvanishing_basis(iu, p, u, U_full);
    
    for l = 1:n_rows
        v = v_vec(l);
        iv = findspan(q, v, n_dimV-1, V_full);
        Nv = nonvanishing_basis(iv, q, v, V_full);
        
        row_idx = (l-1)*n_cols + k;
        
        for jj = 0:q
            for ii = 0:p
                col_idx = (iv-q+jj-1)*n_dimU + (iu-p+ii);
                A(row_idx, col_idx) = Nu(ii+1) * Nv(jj+1);
            end
        end
    end
end

% Resolution of the system:
% If A is ill conditioned, the original problem: min_C ||A*C-Target||
% can be reguralized using Tikonov method:
% min_C ||A*C-Target||+ epsilon ||C||
% => we have to solve the linear system: (A'* A + esp Id)*C = A'* Target
epsilon = 1e-5; % reguralization parameter
ATA = A' * A;
ATb = A' * Target;
C = (ATA + epsilon * speye(size(ATA))) \ ATb;

% compute fitting error: Frobenius norm normalized
Residual = A*C - Target;
err = (norm(Residual, 'fro')^2) / N_tot;
end
%% Example discontinous functions (from Storath and Weimann 2024 sez 3.2)
% WE HAVE DEVELOPED THE STRATEGY EXPLAINED IN THE STORATH AND WEINMANN PAPER
% TO FIND THE DISCONTINUITY SET OF POINTS (J) OF A DISCONTINOUS FUNCTION. 
% AFTER THAT, WE HAVE APPLIED THE KANG'S ALGORITHM ON THE
% NOISY FUNCTION TO FIND THE FINAL ACTIVE KNOTS. 
% WE END UP WITH A SMOOTH SPLINE OF DEGREE 3 OR 5 WHICH INTERPOLATES THE NOISY DATA AND THAT HAS
% DISCONTINUITIES IN J.
% =================================================================

clc; close all; clear;

addpath('../../cvx')
% cvx_setup
addpath('../0. General functions and experiments')

%% parameters
n = 100; a = 0; b = 1; 
p1 = 0.999;               % stiffness parameter p in (0,1) 
noise_level = 0.1;        % standard deviation 
x = linspace(a,b,n)';     % Evaluation points
test_case = 1;            % Choose the function to test

%% =========================================================
% 1. Define test signal
% ==========================================================

switch test_case

    case 1   % -------- g1 --------
        ind1 = (x >= 0.3 & x <= 0.4);
        ind2 = (x >= 0.6);
        f = besselj(1,20*x) + x.*ind1 - x.*ind2; % TRUE SIGNAL

        seed = 10; 
        rng(seed);
        y_data = f + noise_level*randn(n,1);     % NOISY DATA
        gamma = 6;     % hyperparameter selected manually
        tol = 6e-3;    % tollerance defined to get the right plot of the discontinous function
        

    case 2   % -------- g2 (HeaviSine) --------
        f = 4*sin(4*pi*x) - sign(x - 0.3) - sign(0.72 - x); % TRUE SIGNAL
        
        seed = 1; 
        rng(seed);
        y_data = f + noise_level*randn(n,1);                % NOISY DATA
        gamma = 25;    % hyperparameter selected manually
        tol = 6e-3;    % tollerance defined to get the right plot of the discontinous function
end

%% =============================================================
% 2. COMPUTATION OF THE SET OF DISCONTINUITY POINTS J
% ==============================================================

% -------------------------------------------------------------
% 1) Computation of the local energies F*(l,r) by solving the first optimization problem.
% -------------------------------------------------------------

Fstar = inf(n,n);

for l = 1:n
    for r = l+1:n 
        x_loc = x(l:r);
        y_loc = y_data(l:r);
        
        % --- Ar and yr construction---
        [Ar, yr] = build_storath_matrices_local(x_loc, y_loc, p1, noise_level);

        % Least squares problem resolution and computation of the minimum
        % energy as F*(l,r) = min ||Aru-yr||_{2}^{2}
        [Q1,R1] = qr(Ar, 'econ');  
        z = Q1'*yr;
        u = R1\z;
        Fstar(l,r) = sum((Ar * u - yr).^2); 
    end
end

% -------------------------------------------------------------
% 2) Dynamic Programming: Resolution of the second optimization problem.
% Cost vector initializazion: It is inizialized as C = Inf(n,1) because if
% we initialize it as a vector of zeros, then the condition "costo_attuale < 0" 
% would be never satisfied since the energies are all positive or equal
% to zero, so in this case the algorithm would never update the
% compontents' vector.
% -------------------------------------------------------------
C = inf(n,1);
best_l = zeros(n,1);

for r = 1:n
    C(r) = Fstar(1,r);
end

for r = 2:n
    for l = 2:r-2
        % Cost = previous optimal cost + gamma + new energy
        costo_attuale = C(l-1) + gamma + Fstar(l,r);
        if costo_attuale < C(r)
            C(r) = costo_attuale;
            best_l(r) = l;
        end
    end
end

% --------------------------------------------------------
% 3) Backtracking
% --------------------------------------------------------

J_idx = [];
curr = n;
while curr > 0 && best_l(curr) > 0
    split_point = best_l(curr);
    J_idx = [split_point, J_idx];
    curr = split_point - 1;
end

J = (x(J_idx) + x(J_idx-1)) / 2;

fprintf('\nDiscontinuity points: \n')
disp(J)
jump_points = J';

%% REAL FUNCTION PLOT WHICH IS BREAK ON THE DISCONTINUITY POINTS (this is done
% because otherwise matlab connects points in which the function is
% discontinous by default)

f_plot = f; 

for bd = jump_points
    f_plot(abs(x-bd) < tol) = NaN;
end

if test_case == 1
    y_jump_points = -1.1*ones(size(jump_points));
else
    y_jump_points = -7.2*ones(size(jump_points));
end

figure(1); hold on; grid on;
% True signal
plot(x,f_plot,'y','LineWidth',2)
% Noisy data
plot(x, y_data, 'm.'); 
% Discontinuity points
plot(jump_points,y_jump_points,'rv','MarkerFaceColor','r');
title('True discontinuous signal')
legend('True signal', 'Noisy data', 'Discontinuity points')
xlabel('x')

%% =========================================================
% 3. KANG ALGORITHM APPLICATION
% =========================================================

% ------------------------------------------------------
% 1) Set up parameters and Initial inner knot vector definitions
% ------------------------------------------------------
% Spline degree
p = input('Insert the degree of the spline p (3 or 5): \n');

if p ~= 3 && p ~= 5
    error('Insert the correct degree: p=3 or p=5')
end

N = 100;  
n_init_knots = 40;    % Initial number of knots
inner_knots = linspace(a,b,n_init_knots); % Initial knots points
inner_knots = inner_knots(2:end-1);       % Remove the end points, consider only the inner knots

% impose jumps points multiplicity = p+1, so that the reconstructed
% B-spline has discontinuities in J.
jump_knots = repelem(jump_points, p+1);

% Initial knot vector (CLAMPED)
U_init = [a*ones(1,p+1), inner_knots, jump_knots, b*ones(1,p+1)];
U_init = sort(U_init);

% RECALL: n+1 = m-p => n = m-p-1 where U = (u0=...=up, ..., um-p=...=um) (Recall in matlab you have to shift
% indices)
n_dim = length(U_init) - p - 1;
n_kang = n_dim - 1;

% --------------------------------------------------------------
% 2.1). Construction of A: A c = y
% The matrix A has dimension N x n_dim => each row k contains the values of
% the basis fcts Ni,p evaluated in x(k).
% --------------------------------------------------------------

% Pre allocation for A (sparse format)
A_kang = sparse(N,n_dim);
% Construct A row-wise (for each element x(k))
for k = 1:N
    ik = findspan(p,x(k),n_kang,U_init);  % Find the knot span of U_init where x(k) lives
    A_kang(k,ik-p:ik) = nonvanishing_basis(ik,p,x(k),U_init);   % Find the non vanishing basis fcts (evaluated in x(k))
    % This fct returns a p+1 vector containing the p+1 non vanishing basis in
    % the knot span of x(k), and save these values of the fcts in the k th row of A and corresponding columns (i_k-p ... i_k)
end

% -----------------------------------------------------------------
% 2.2). Construction of the jump matrix D.
% We need to compute in each inner node the jump of the pth order
% derivative of the curve, but the pth order derivative of the curve is
% given by the lin. comb. of c times the pth derivative of the basis fcts
% Ni,p.
% -----------------------------------------------------------------

delta_num = 1e-4;
num_jumps = length(inner_knots);

% Pre allocate D
% Each row of D corresponds to a different internal node, while the columns
% correspond to the basis fcts
D = zeros(num_jumps,n_dim);
% Loop over each inner node for computing the jump
for j = 1:num_jumps
    u = inner_knots(j);
    % Right limit (u + delta)
    ikR = findspan(p,u+delta_num,n_kang,U_init);
    DR  = derivative_basis(ikR,p,u+delta_num,U_init,p);
    % The fct return a matrix containing in each row the derivative of
    % order k of the basis fcts non vanishing at the given knot span
    
    % Left limit (u - delta)
    ikL = findspan(p,u-delta_num,n_kang,U_init);
    DL  = derivative_basis(ikL,p,u-delta_num,U_init,p);
    % Insert the values corresponding to the derivatives of order p+1 
    % on the jth row and corresponding columns: 
    D(j,ikR-p:ikR) = D(j,ikR-p:ikR) + DR(p+1,:);
    % Subtract to obtain the jump: 
    D(j,ikL-p:ikL) = D(j,ikL-p:ikL) - DL(p+1,:);
end

% Normalize the constructed matrix D in order to scale it, since for p = 5
% it will have components of high degree
D = D./max(abs(D(:)));

% --------------------------------------------------------
% 3) STEP 1: Sparse optimization + selection active knots
% ---------------------------------------------------------
if p == 3
    lambda_kang = 1e-6; % parameter chosen sperimentally used in the definition of the convex optimization problem 
else
    lambda_kang = 1e-3; 
end

% Start a block to use the CVX package 
% compute c as the solution of the convex optimization problem
cvx_begin quiet
    variable c(n_dim)     % n_dim length of the vector of the variables c
    
    % Minimize the norm of Ac-f (in this case this is a vector so we
    % use the 2 norm) + lambda * norm 1 of Dc
    minimize( norm(A_kang*c - y_data,2) + lambda_kang*norm(D*c,1) )
cvx_end

% Find active knots
% thereshold selection
if p == 3
    thereshold = 1e-3;
else
    thereshold = 1e-6;
end

jump_measure = abs(D*c);
idx_active = find(jump_measure > thereshold);   % Indices corresponding to the active knots: the knots where the jump is non zero
U_active = inner_knots(idx_active);             % Extract the active knots from the inner knots

% --------------------------------------------------------
% 4) STEP 2 – Adjustment & pruning
% --------------------------------------------------------
% Tolerance for the adjustment of the knots
tol_adj = 1e-4;
% Parameter for constructing clusters chosen sperimentally
delta_cluster = 2*(inner_knots(2) - inner_knots(1));
active_knots_adj = AdjustKnotsPosition( delta_cluster, tol_adj, U_active, p, x, y_data, a, b);  

fprintf('\nFinal inner knots: \n')
disp(active_knots_adj);
fprintf('Number of final inner knots: %d\n', length(active_knots_adj));

%% =========================================================
% 5) Final currupted spline reconstruction using the final adjusted knots
% ==========================================================
% Final knot vector
U_final = [a*ones(1,p+1), active_knots_adj, jump_knots, b*ones(1,p+1)];
U_final = sort(U_final);
% Compute the dimension to construct A:
nF = length(U_final)-p-1;

A_final = sparse(N,nF);
for k = 1:N
    ik = findspan(p,x(k),nF-1,U_final);
    A_final(k,ik-p:ik) = nonvanishing_basis(ik,p,x(k),U_final);
end
% Solve the final Least square pb:
c_final = A_final\y_data;

curve_final = zeros(N,1);
for k = 1:N
    ik = findspan(p,x(k),nF-1,U_final);
    B  = nonvanishing_basis(ik,p,x(k),U_final);
    curve_final(k) = B*c_final(ik-p:ik);
end

% Break the reconstructed curve at the discontinuity points
curve_plot = curve_final;
for bb = jump_points
    curve_plot(abs(x-bb)<tol) = NaN;
end

%% Final plot
% Comparison between the true signal and the reconstructed cubic spline approximation.
if test_case == 1
    y_knots_adj = -1.1*ones(size(active_knots_adj));
    y_jump_points = -1.1*ones(size(jump_points));
else
    y_knots_adj = -7.2*ones(size(active_knots_adj));
    y_jump_points = -7.2*ones(size(jump_points));
end

figure(2); hold on; grid on;
% True signal
plot(x,f_plot,'r','LineWidth',2)
% Reconstructed corrupted signal
plot(x,curve_plot,'b','LineWidth',2)
% Adjusted final knots
plot(active_knots_adj, y_knots_adj,'g^','MarkerFaceColor','g')
% Discontinuity points
plot(jump_points, y_jump_points,'rv','MarkerFaceColor','r');
% Noisy data
plot(x, y_data, 'm.');
title(sprintf('Final B-spline reconstruction (p = %d)',p))
legend('True signal','Reconstructed currupted spline','Adaptive knots', 'Discontinuous points', 'Noisy data')

%% MSE and ME
% Residual:
res = A_final*c_final - f;

% MSE: mean square error 
MSE = sum(res.^2) / N;

% ME: maximum error
ME = max(abs(res));

fprintf('MSE: %e\n', MSE);
fprintf('ME:  %e\n', ME);
%% Example Parametric Curve - Cusp adding some Noise
% From Luo et al. (2021) - section 4.3

clc; 
close all; 
clear;

addpath('../../cvx')
% cvx_setup
addpath('../0. General functions and experiments')

%% Construction of the Parametric Open curve

a=0;
b=1;

N = 1001;  % Number of data points                           
t = linspace(a,b,N)';  %parameter t

x_true = abs(4*t - 2) + 1;
y_true = (4*t - 2).^3 + t;
P_true = [x_true, y_true];

% Add of Noise
rng(42)
delta_noise = 0.08; 
% For each point (x,y) we add a random number sampled uniformly in [-0.02 0.02]
P_noisy = P_true + delta_noise * (2 * rand(N, 2) - 1);

%% Setup parameters

% Spline degree
p = input('Insert the degree of the spline p (3 or 5): \n');

if p ~= 3 && p ~= 5
    error('Insert the correct degree: p=3 or p=5')
end                           

% lambda regularization parameter
if p == 3
    lambda = 1e-5;
    
else
    lambda = 1e-10;
end

n_init_knots = 101;                       % Initial number of knots
inner_knots = linspace(a,b,n_init_knots); % Initial knot points
inner_knots = inner_knots(2:end-1);       % Remove the end points, consider only the inner knots

% Initial knot vector (CLAMPED);
U_init = [a*ones(1,p+1), inner_knots, b*ones(1,p+1)];

% RECALL: n+1 = m-p => n = m-p-1 where U = (u0=...=up, ..., um-p=...=um) (Recall in matlab you have to shift
% indices) 
n_dim = length(U_init) - p - 1; %n+1
n = n_dim - 1;

%% Step 1
%% Construction of A:
% The matrix A has dimension Nxn_dim => each row k contains the values of
% the basis fcts Ni,p evaluated at t in [0,1]

% Pre allocation for A (sparse format)
A = sparse(N,n_dim);

% Construct A row-wise (for each element t(k))
for k = 1:N
    ik = findspan(p,t(k),n,U_init); % Find the knot span of U_init where t(k) lives

    basis_vals = nonvanishing_basis(ik, p, t(k), U_init);     % Find the non vanishing basis fcts (evaluated in x_norm(k))
    % This fct returns a p+1 vector containing the p+1 non vanishing basis in
    % the knot span of x_norm(k)

    A(k,ik-p:ik) = basis_vals;  % Save these values of the fcts in the k-th row of A and corresponding columns (i_k-p ... i_k)
end

%% Construction of D the jump matrix
% We need to compute in each inner node the jump of the pth order
% derivative of the curve, but the pth order derivative of the curve is
% given by the lin. comb. of c times the pth derivative of the basis fcts
% Ni,p

% Consider the unique knots (= inner knots) for the jumps:
unique_inner_knots = inner_knots; 
num_jumps = length(unique_inner_knots);

% Pre allocate D
% Each row of D corresponds to a different internal node, while the columns
% correspond to the basis fcts
D = zeros(num_jumps,n_dim);

% Choose a small delta for computing the jump (right and left derivative at
% a point u)
delta = 1e-6;

% Loop over each inner node for computing the jump
for j = 1:num_jumps
    u_node = unique_inner_knots(j);
    
    % Right limit (u + delta)
    ik_R = findspan(p, u_node + delta, n, U_init);
    DN_R  = derivative_basis(ik_R, p, u_node + delta, U_init, p);
    % The fct return a matrix containing in each row the derivative of
    % order k of the basis fcts non vanishing at the given knot span

    % Insert the values corresponding to the derivatives of order p+1 
    % on the jth row and corresponding columns: 
    % (Recall D at the begininning is a matrix of zeros)
    D(j,ik_R-p:ik_R) = D(j,ik_R-p:ik_R) + DN_R(p+1,:);
    
    % Left limit (u - delta)
    ik_L = findspan(p, u_node - delta, n, U_init);
    DN_L  = derivative_basis(ik_L, p, u_node-delta, U_init, p);
    % Subtract to obtain the jump: 
    D(j,ik_L-p:ik_L) = D(j,ik_L-p:ik_L) - DN_L(p+1,:);
end

% Start a block to use the CVX package

% cvx_begin quiet           % (if you don't want everything printed)
cvx_begin
    variable C(n_dim,2)
    minimize( norm(A*C - P_noisy,'fro') + lambda*sum(norms(D*C,Inf,2)))
cvx_end

%% Find the active knots
jump_norm = sqrt(sum((D*C).^2,2));

idx_active = find(jump_norm > 10);    %Indices corresponding to the active knots: the knots where the jump is non zero
active_knots = inner_knots(idx_active); % Extract the active knots from the inner knots

U_active = [a*ones(1, p+1), active_knots, b*ones(1, p+1)];        % Construct the clamped knot vector with the new inner knots

fprintf('Knots after Step 1: %d\n', length(active_knots));
disp(active_knots);

%% Step 2:

% Define delta cluster 
delta_cluster = inner_knots(2)-inner_knots(1);

% Tolerance for the adjustment of the knots:
tol_adj = 1e-3;

% Adjust the final internal knot:
active_knots_adj = AdjustKnotsPosition(delta_cluster, tol_adj, active_knots, p, t, P_noisy, a, b);

fprintf('Knots after Step 2: %d\n', length(active_knots_adj));
disp(active_knots_adj);

%% Final spline reconstruction:

U_final_adj = [a*ones(1,p+1), active_knots_adj, b*ones(1,p+1)];

% Compute the dimension to construct A:
n_dim_final = length(U_final_adj) - p - 1;
n_final = n_dim_final - 1;

% Construct the matrix A:
A_final = sparse(N, n_dim_final);
for k = 1:N
    ik = findspan(p, t(k), n_final, U_final_adj);
    basis_vals = nonvanishing_basis(ik, p, t(k), U_final_adj);
    A_final(k,(ik-p):ik) = basis_vals;
end

% Solve the final Least square pb:

C_final = A_final \ P_noisy;

% Final plot:
tt = linspace(a,b,N);
curve_final = zeros(length(tt),2);
for k = 1:length(tt)
    ik = findspan(p, tt(k), n_final, U_final_adj);
    B  = nonvanishing_basis(ik, p, tt(k), U_final_adj);
    curve_final(k,:) = B*C_final(ik-p:ik,:);
end

x_knot = abs(4*active_knots_adj - 2) + 1;
y_knot = (4*active_knots_adj - 2).^3 + active_knots_adj;
P_knot = [x_knot', y_knot'];

figure(1); 
hold on; 
grid on;
plot(P_noisy(:,1),P_noisy(:,2),'ro');
plot(curve_final(:,1),curve_final(:,2),'b-','LineWidth',2);
plot(P_true(:,1),P_true(:,2),'k--','LineWidth',1.5);
plot(P_knot(:,1),P_knot(:,2),'g^','LineWidth',1.5)
legend('Given data noised','B-spline approximation','Original data','Final knots');
title(sprintf('Final B-spline reconstruction (p = %d)',p));
xlim([0.8, 3.5])

%% MSE and ME
% Residual:
res = A_final * C_final - P_noisy;
error = vecnorm(res,2,2);  % Returns a vector of the norm 2 along the 2nd component (column),
                           % so for each row computes the 2-norm: sqrt(resx^2 + resy^2)

ME = max(error);
MSE = (norm(res,'fro')^2) / N;

fprintf('ME:  %e\n', ME);
fprintf('MSE: %e\n', MSE);
%% Example: Extension of Kang with noisy data and DE algorithm
% From Luo et al. (2021) - section: 4.3
clc; 
clear; 
close all;

addpath('../../cvx')
% cvx_setup
addpath('../0. General functions and experiments')
addpath('differentialevolution 2019-11-13') 

%% Construction of the Parametric Open curve:

% Flag variable for adding or not noise:
test_case_noise = true; % false = clean data, true = noisy data


% Domain end points:
a = 0; 
b = 1;

N = 51;          % Data points 
t = linspace(a, b, N)';

x_true = abs(4*t - 2) + 1;
y_true = (4*t - 2).^3 + t;
P_true = [x_true, y_true];

% Generation of the data points (Target curve)
P_target = [abs(4*t - 2) + 1, (4*t - 2).^3 + t];

% Choose if we need to add noise or not:
rng(42)     % Fix the seed
if test_case_noise
    % Add noise: 
    delta_noise = 0.02;
    % For each point (x,y) we add a random number sampled uniformly in [-0.02 0.02]
    P_target = P_target + delta_noise * (2*rand(N,2) - 1);
end

%% Setup parameters:
% Reconstructed spline degree:
p = input('Insert the degree of the spline p (3 or 5): \n');

if p ~= 3 && p ~= 5
    error('Insert the correct degree: p=3 or p=5')
end  

% lambda = regularization parameter
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
cvx_begin quiet
    variable C(n_dim, 2)
    minimize( norm(A*C - P_target, 'fro') + lambda*sum(norms(D*C, Inf, 2)) )
cvx_end

%% CONTROLL
% % --- NUOVA LOGICA DI SELEZIONE: I 3 PICCHI PIÙ SIGNIFICATIVI ---
jump_norm = sqrt(sum((D*C).^2, 2));

% [~, sorted_idx] = sort(jump_norm, 'descend');
% num_nodi_da_cercare = 3; % Il paper suggerisce che bastano pochi nodi per la cuspide
% idx_active = sorted_idx(1:num_nodi_da_cercare);
% U_active = sort(inner_knots(idx_active));

if p == 3
    eta = 5e-3;                 % Relative threshold
else 
    eta = 2e-1;
end

% Algorithm 3.1 LUO 2021
% Find the local maxima:
K0 = []; 
j = 1;
for i = 2:num_jumps-1
    jump_lmax = max([jump_norm(i-1), jump_norm(i), jump_norm(i+1)]);
    if jump_lmax == jump_norm(i)     % Local maximum
        K0(j) = i; 
        j = j + 1;
    end
end

% filtro per soglia relativa
K1 = find(jump_norm(K0) > eta * max(jump_norm));  % Select the indices inside K0 whose correspondent jump is above the rel threshold
idx_active = K0(K1);                              % True indices on inner_knots

U_active = sort(inner_knots(idx_active));            % Active knots


fprintf('Knots after Step 1: %d\n', length(U_active));
disp(U_active);


%% Step 2:
% Knots optimization by DE algorithm
% Construct a structure to give in input to the fittness function:
objFctSettings.P = P_target;        % Save the data points
objFctSettings.deg = p; % Save the spline degree
objFctSettings.t = t;           % Save the parameter vector t associated to the data points
objFctSettings.a = a;           % Save the end points of the parametric domain
objFctSettings.b = b;

% Construct a structure with the default parameters of the DE library:
DEParams = getdefaultparams;

DEParams.NP = 15 * length(U_active);    % Fix the population number: Each element of the 
% population is a vector of length = length(U_active) and we choose L = 15
% (as suggested in the paper in [10,20])

DEParams.maxiter = 100;       % Fix the max number of iterations
DEParams.displayResults = 0;   % We don't want to display the results
DEParams.saveHistory = 0;       % Don't save in the memory the complete history (in this way we save memory)

% Construct an array of dimension Kx4 (with K = lentgh(U_active)
% Each row  describes an optimized parameter (knot)
paramDefCell = cell(length(U_active), 4);
% Loop over the active knots (that DE will move):
for i = 1:length(U_active)
    paramDefCell{i,1} = sprintf('u%d', i);      % Symbolic name of the parameter u1, u2...
    paramDefCell{i,2} = [a + 0.05, b - 0.05];   % Bounds for the search space
    paramDefCell{i,3} = 0;      % Flag variable to define the type of variable (0=continuous)
    paramDefCell{i,4} = U_active(i);        % Initial guess for the parameter (knot)
end

% Call the DE function 
% It returns the "best element" found = optimized position of the knots

[bestmem, ~] = differentialevolution(DEParams, paramDefCell, @fitness_function, objFctSettings);
U_adjusted = sort(bestmem);     % Order the values given by DE

fprintf('Knots after Step 2: %d\n', length(U_adjusted));
disp(U_adjusted');

%% Final spline reconstruction:
U_final = [a*ones(1, p+1), U_adjusted(:)', b*ones(1, p+1)];

% Compute the dimension to construct A
n_dim_final = length(U_final) - p - 1;
n_final = n_dim_final - 1;

% Construct the matrix A:
A_final = sparse(N, n_dim_final);
for k = 1:N
    ik = findspan(p, t(k), n_final, U_final);
    basis_vals = nonvanishing_basis(ik, p, t(k), U_final);
    A_final(k,(ik-p):ik) = basis_vals;
end

C_final = A_final \ P_target;

curve_final = A_final * C_final;

x_knot = abs(4*U_adjusted - 2) + 1;
y_knot = (4*U_adjusted - 2).^3 + U_adjusted;
P_knot = [x_knot(:), y_knot(:)];

figure(1); 
hold on; 
grid on;
plot(P_target(:,1),P_target(:,2),'ro');
plot(curve_final(:,1),curve_final(:,2),'b-','LineWidth',2);
plot(P_true(:,1),P_true(:,2),'k--','LineWidth',1.5);
plot(P_knot(:,1),P_knot(:,2),'g^','LineWidth',1.5)
legend('Given data noised','B-spline approximation','Original data','Final knots');
title(sprintf('Final B-spline reconstruction (p = %d)',p));
xlim([0.8, 3.5])


%% MSE and ME
% Residual:
res = A_final * C_final - P_target;
error = vecnorm(res,2,2);  % Returns a vector of the norm 2 along the 2nd component (column),
                           % so for each row computes the 2-norm: sqrt(resx^2 + resy^2)

ME = max(error);
MSE = (norm(res,'fro')^2) / N;

fprintf('ME:  %e\n', ME);
fprintf('MSE: %e\n', MSE);


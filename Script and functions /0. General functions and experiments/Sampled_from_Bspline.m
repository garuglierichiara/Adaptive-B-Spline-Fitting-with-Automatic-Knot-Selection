%% Example: Sampled from a B-spline function (from Kang 2015)
clc; 
close all; 
clear;

addpath('../../cvx')
% cvx_setup

%% Construction of the Original Function

%spline degree fixed for construction of the Original Function
p_true = 3;

% Parameters
a = 0; 
b = 1; 
N = 1001;         %Number of Data points

% Given internal knots 
U_true_int = [0.0439, 0.0653, 0.2293, 0.2367, 0.4821, 0.4907, 0.5408, 0.5408, 0.6209, 0.7051, 0.9407];

% Initial knot vector (CLAMPED) U = (u0=...=up, ..., um-p=...=um)
U_true = [a*ones(1,p_true+1), U_true_int, b*ones(1,p_true+1)];

% Number of Control points: follow the rule m - p = n + 1 
% where m+1 = length(U_init), p = degree
n_dim_true = length(U_true) - p_true - 1; %n+1
Ptrue = [-15;10;37;47;-52;80;43;76;10;3;69;26;97;80;78];
n_true = length(Ptrue)-1;

x = linspace(a, b, N)';
f = zeros(size(x));
for i = 1:length(x)
    f(i) = bspline_evaluation(p_true, x(i), n_true, U_true, Ptrue); 
end

figure(1)
hold on;
y = -18 * ones(size(U_true_int)); 
plot(U_true_int, y, 'b^', 'MarkerSize', 8,'MarkerFaceColor', 'b');
plot(x,f,'r')
title('Original function built with spline of degree 3')
xlim([0, 1]);
ylim([-20,100]);
legend('Initial knots', 'Original function');
grid on;

%% Setup parameters

%epsilon = 2e-5;    % Tolerance
lambda = 0.5;       % To choose sperimentally

% Spline degree to approximate the original function
p = input('Insert the degree of the spline p (3 or 5): \n');

if p ~= 3 && p ~= 5
    error('Insert the correct degree: p=3 or p=5')
end

n_init_knots = 501;       % Initial number of knots
inner_knots = linspace(a, b, n_init_knots); % Initial knot points
inner_knots = inner_knots(2:end-1); %togliamo gli estremi, consideriamo solo i nodi interni

% Initial knot vector (CLAMPED);
U_init = [a*ones(1, p+1), inner_knots, b*ones(1, p+1)];

% RECALL: n+1 = m-p => n = m-p-1 where U = (u0=...=up, ..., um-p=...=um) (Recall in matlab you have to shift
% indices) 
n_dim = length(U_init)- p - 1;  % = n+1
n = n_dim-1;

%% Step 1
%% Construction of A: A c = y_norm
% The matrix A has dimension Nxn_dim => each row k contains the values of
% the basis fcts Ni,p evaluated in x_norm(k)

% Pre allocation for A (sparse format)
A = sparse(N, n_dim); 

% Construct A row-wise (for each element x_norm(k))
for k = 1:N
    i_k = findspan(p, x(k), n, U_init);    % Find the knot span of U_init where x_norm(k) lives

    basis_vals = nonvanishing_basis(i_k, p, x(k), U_init);     % Find the non vanishing basis fcts (evaluated in x_norm(k))
    % This fct returns a p+1 vector containing the p+1 non vanishing basis in
    % the knot span of x_norm(k)

    A(k, (i_k-p):i_k) = basis_vals;     % Save these values of the fcts in the k th row of A and corresponding columns (i_k-p ... i_k)
end


%% Construction of D the jump matrix
% We need to compute in each inner node the jump of the pth order
% derivative of the curve, but the pth order derivatice of the curve is
% given by the lin. comb. of c times the pth derivative of the basis fcts
% Ni,p

% Consider the unique knots (= inner knots) for the jumps:
unique_inner_knots = inner_knots; 
num_jumps = length(unique_inner_knots);

% Pre allocate D
% Each row of D corresponds to a different internal node, while the columns
% correspond to the basis fcts
D = zeros(num_jumps, n_dim);

% Choose a small delta for computing the jump (right and left derivative at
% a point u)
delta = 1e-6; 

% Loop over each inner node for computing the jump
for j = 1:num_jumps
    u_node = unique_inner_knots(j);

    % Right limit (u + delta)
    ik_R = findspan(p, u_node + delta, n, U_init);
    DN_R = derivative_basis(ik_R, p, u_node + delta, U_init, p);       
    D(j, (ik_R-p):ik_R) = D(j, (ik_R-p):ik_R) + DN_R(p + 1, :);

    % Left limit (u - delta)
    ik_L = findspan(p, u_node - delta, n, U_init);
    DN_L = derivative_basis(ik_L, p, u_node - delta, U_init, p);
    % Subtract to obtain the jump: 
    D(j, (ik_L-p):ik_L) = D(j, (ik_L-p):ik_L) - DN_L(p + 1, :);
end

% Normalize D: to make comparable the norm of A and D,
% Otherwhise we won't get a good balance between the minimization problem
% of the fidelity term vs minimization of the regularization term 
D = D ./ max(abs(D(:)));

% Start a block to use the CVX package

%cvx_begin quiet            % (if you don't want everything printed)
cvx_begin
    variable c(n_dim) % n_dim length of the vector of the variables c

    % Minimize the norm of Ac-y_norm (in this case this is a vector so we
    % use the 2 norm) + lambda * norm 1 of Dc
    minimize( norm(A*c - f, 2) + lambda * norm(D*c, 1))
cvx_end

%% Prepare the results for the plot

% Generate 1000 points to plot the curve
xx = linspace(a, b, 1000);
curve = zeros(size(xx));

for k = 1:length(xx)
    % Use findspan and the coefficients "c" computed by cvx
    idx = findspan(p, xx(k), n, U_init);
    basis_fcts = nonvanishing_basis(idx, p, xx(k), U_init);

    % curve = linear combination of basis_fcts*coefficients c
    curve(k) = basis_fcts * c(idx-p:idx);
end

%% Find the active knots
jump = D*c;

idx_active = find(abs(jump)>3*1e-5); % Indices corresponding to the active knots: the knots where the jump is non zero

active_knots = inner_knots(idx_active); % Extract the active knots from the inner knots

U_active = [a*ones(1, p+1), active_knots, b*ones(1, p+1)];        % Construct the clamped knot vector with the new inner knots

fprintf('\nActive knots: \n'); 
disp(active_knots);

% Plot
figure(2); 
hold on;
grid on;
plot(x, f, 'r', 'LineWidth', 2.5); 
plot(xx, curve, 'y--', 'LineWidth', 1.5);
y_knots = -18 * ones(size(active_knots)); 
plot(active_knots, y_knots, 'b^', 'MarkerSize', 6,'MarkerFaceColor', 'b');
title(sprintf('Initial B-spline reconstruction and active knots (p = %d)', p));
ylim([-20, 100]); 
xlim([0, 1]);
legend('Original data', 'Reconstructed spline', 'Active knots');

%% Step 2:

% Define delta cluster
%delta_cluster = 0.008;
delta_cluster = 2*(inner_knots(2) - inner_knots(1));

% Tolerance for the adjustment of the knots:
tol_adj = 1e-5;

% Adjust the final internal knot:
active_knots_adj = AdjustKnotsPosition(delta_cluster, tol_adj, active_knots, p, x, f, a, b);

fprintf('\nFinal inner knots: \n'); 
disp(active_knots_adj);

%% Step 3:
U_final_adj = [a*ones(1, p+1), active_knots_adj, b*ones(1, p+1)];

% Compute the dimension to construct A:
n_dim_final = length(U_final_adj) - p - 1;
n_final = n_dim_final - 1;

% Construct the matrix A:
A_final = sparse(N, n_dim_final);
for k = 1:N
    ik = findspan(p, x(k), n_final, U_final_adj);
    basis_vals = nonvanishing_basis(ik, p, x(k), U_final_adj);
    A_final(k, (ik-p):ik) = basis_vals;
end

% Solve the final Least square pb:
c_final = A_final \ f;

% Final plot:
curve_final = zeros(size(xx));
for k = 1:length(xx)
    idx = findspan(p, xx(k), n_final, U_final_adj);
    b_f = nonvanishing_basis(idx, p, xx(k), U_final_adj);
    curve_final(k) = b_f * c_final(idx-p:idx);
end

% Plot
figure(3); 
hold on;
plot(x, f, 'ro', 'MarkerFaceColor','red');
plot(xx, curve_final, 'y-', 'LineWidth', 2);
y_knots_adj = -18* ones(size(active_knots_adj)); 
plot(active_knots_adj, y_knots_adj, 'b^', 'MarkerSize', 8,'MarkerFaceColor', 'b');
title(sprintf('Final B-spline reconstruction (p = %d)',p));
legend('Original data','Reconstructed spline', 'Final knots');
xlim([0, 1]);
ylim([-20,100]);
grid on;

%% MSE and ME
% Residual:
res = A_final*c_final - f;

% MSE: mean square error 
MSE = sum(res.^2) / N;

% ME: maximum error
ME = max(abs(res));

fprintf('MSE: %e\n', MSE);
fprintf('ME:  %e\n', ME);

%% Example: Titanium heat data 
% From Kang et al. (2015)
clc;
close all;
clear;

addpath('../../cvx')
% cvx_setup

% Extract Titanium heat data (From the Curve Fitting Toolbox)
[x, y] = titanium;

% Spline degree
p = input('Insert the degree of the spline p (3 or 5): \n');

if p ~= 3 && p ~= 5
    error('Insert the correct degree: p=3 or p=5')
end


% Plot of the original data points:
figure(1);
plot(x,y,'ro', MarkerFaceColor='red');
axis([580, 1090, 0.4, 2.3]);
title('Titanium heat data');
legend('Original data');
grid on;

%% Normalization of the data:
% New interval:
a = 0;
b = 75;

x_min = min(x);
x_max = max(x);

% Normalize the dataset in [0, 75]
% Trasliamo e riscaliamo la curva per definirla su [a,b] 
x_norm = a + (x - x_min) * (b - a) / (x_max - x_min);
y_norm = y(:);

% Plot of the normalized data:
figure(2);
plot(x_norm,y_norm,'ro', MarkerFaceColor='red');
title('Titanium heat data normalized');
axis([-5, 80, 0.4, 2.3]);
legend('Original data');
grid on;


%% Setup parameters:

N = length(x_norm);       % Number of data points
n_init_knots = 101;       % Initial number of knots             

if p==3
    lambda = 0.5;
else                      % Regularization parameter
    lambda = 10;
end

inner_knots = linspace(a, b, n_init_knots); % Initial knot points
inner_knots = inner_knots(2:end-1);         % Remove the end points, consider only the inner knots 

% Initial knot vector (CLAMPED);
U_init = [a*ones(1, p+1), inner_knots, b*ones(1, p+1)];

% RECALL: n+1 = m-p => n = m-p-1 where U = (u0=...=up, ..., um-p=...=um) (Recall in matlab you have to shift
% indices) 
n_dim = length(U_init)- p - 1;  % = n+1

n = n_dim-1;

%% Step 1:
%% Construction of A: A c = y_norm
% The matrix A has dimension Nxn_dim => each row k contains the values of
% the basis fcts Ni,p evaluated in x_norm(k)

% Pre allocation for A (sparse format)
A = sparse(N, n_dim); 

% Construct A row-wise (for each element x_norm(k))
for k = 1:N
    i_k = findspan(p, x_norm(k), n, U_init);    % Find the knot span of U_init where x_norm(k) lives
    
    basis_vals = nonvanishing_basis(i_k, p, x_norm(k), U_init);     % Find the non vanishing basis fcts (evaluated in x_norm(k))
    % This fct returns a p+1 vector containing the p+1 non vanishing basis in
    % the knot span of x_norm(k)

    A(k, (i_k-p):i_k) = basis_vals;     % Save these values of the fcts in the k th row of A and corresponding columns (i_k-p ... i_k)
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
    % The fct return a matrix containing in each row the derivative of
    % order k of the basis fcts non vanishing at the given knot span
    
    % Insert the values corresponding to the derivatives of order p+1 
    % on the jth row and corresponding columns: 
    % (Recall D at the begininning is a matrix of zeros)
    D(j, (ik_R-p):ik_R) = D(j, (ik_R-p):ik_R) + DN_R(p + 1, :);
    
    % Left limit (u - delta)
    ik_L = findspan(p, u_node - delta, n, U_init);
    DN_L = derivative_basis(ik_L, p, u_node - delta, U_init, p);
    % Subtract to obtain the jump: 
    D(j, (ik_L-p):ik_L) = D(j, (ik_L-p):ik_L) - DN_L(p + 1, :);
end

% Normalization of D in case of degree 5
if p==5
    D = D ./ max(abs(D(:)));
end

% Start a block to use the CVX package

% cvx_begin quiet           % (if you don't want everything printed)
cvx_begin
    variable c(n_dim)           % n_dim length of the vector of the variables c
    
    % Minimize the norm of Ac-y_norm (in this case this is a vector so we
    % use the 2 norm) + lambda * norm 1 of Dc
    minimize( norm(A*c - y_norm, 2) + lambda * norm(D*c, 1))     % Objective function
cvx_end

%% Prepare the results for the plot
 
% Generate 1000 points to plot the curve
xx = linspace(a, b, 1000);
curve1 = zeros(size(xx));
    
for k = 1:length(xx)
    % Use findspan and the coefficients "c" computed by cvx
    idx = findspan(p, xx(k), n, U_init);
    basis_fcts = nonvanishing_basis(idx, p, xx(k), U_init);
    
    % curve = linear combination of basis_fcts*coefficients c
    curve1(k) = basis_fcts * c(idx-p:idx);
end

%% Find the active knots
jump = D*c;

if p==3
    idx_active = find(abs(jump)>1e-3);      % Indices corresponding to the active knots: the knots where the jump is non zero
else
    idx_active = find(abs(jump)>1e-4);      % Indices corresponding to the active knots: the knots where the jump is non zero
end
active_knots = inner_knots(idx_active); % Extract the active knots from the inner knots

U_active = [a*ones(1, p+1), active_knots, b*ones(1, p+1)];        % Construct the clamped knot vector with the new inner knots

fprintf('\nActive knots: \n'); 
disp(active_knots);

% Plot
figure(3); 
hold on;
grid on;
plot(x_norm, y_norm, 'ro', 'MarkerFaceColor','red'); 
plot(xx, curve1, 'b-', 'LineWidth', 2);
y_knots = 0.5 * ones(size(active_knots)); 
% Plot the active knots as blue triangles
plot(active_knots, y_knots, 'b^', 'MarkerSize', 8,'MarkerFaceColor', 'b');
legend('Original data points','Reconstructed spline', 'Active knots');
title(sprintf('Initial B-spline reconstruction and active knots (p = %d)', p));
axis([-5, 80, 0.4, 2.3]);

%% Step 2:

% Define delta cluster as the "length of two adjacent inner knots in the initial knot vector U": 
delta_cluster = inner_knots(2) - inner_knots(1); 

% Tolerance for the adjustment of the knots:
tol_adj = 1e-3;      

% Adjust the final internal knot:
active_knots_adj = AdjustKnotsPosition(delta_cluster, tol_adj, active_knots, p, x_norm, y_norm, a, b);

fprintf('\nFinal inner knots: \n'); 
disp(active_knots_adj);

%% Final spline reconstruction:
U_final_adj = [a*ones(1, p+1), active_knots_adj, b*ones(1, p+1)];

% Compute the dimension to construct A:
n_dim_final = length(U_final_adj) - p - 1;
n_final = n_dim_final - 1;

% Construct the matrix A:
A_final = sparse(N, n_dim_final);
for k = 1:N
    ik = findspan(p, x_norm(k), n_final, U_final_adj);
    basis_vals = nonvanishing_basis(ik, p, x_norm(k), U_final_adj);
    A_final(k, (ik-p):ik) = basis_vals;
end

% Solve the final Least square pb:
c_final = A_final \ y_norm;

% Final plot:
curve_final = zeros(size(xx));
for k = 1:length(xx)
    idx = findspan(p, xx(k), n_final, U_final_adj);
    b_f = nonvanishing_basis(idx, p, xx(k), U_final_adj);
    curve_final(k) = b_f * c_final(idx-p:idx);
end

figure(4);
hold on;
plot(x_norm, y_norm, 'ro', 'MarkerFaceColor','red');    % Original data points
plot(xx, curve_final, 'b-', 'LineWidth', 2);            % Approximated spline
y_knots_adj = 0.5* ones(size(active_knots_adj)); 
plot(active_knots_adj, y_knots_adj, 'b^', 'MarkerSize', 8,'MarkerFaceColor', 'b');
legend('Original data points','Reconstructed spline','Final knots');
title(sprintf('Final B-spline reconstruction (p = %d)',p));
axis([-5, 80, 0.4, 2.3]);
grid on;

%% MSE and ME
% Residual:
res = A_final*c_final - y_norm;

% MSE: mean square error 
MSE = sum(res.^2) / N;

% ME: maximum error
ME = max(abs(res));

fprintf('MSE: %e\n', MSE);
fprintf('ME:  %e\n', ME);

% Residual error defined as in the paper
w = ones(N,1);
w(1) = 1/2; w(end) = 1/2;

res_err = sqrt(sum(w.*(res.^2))/(N-1));

fprintf('Residual error:  %e\n\n', res_err);
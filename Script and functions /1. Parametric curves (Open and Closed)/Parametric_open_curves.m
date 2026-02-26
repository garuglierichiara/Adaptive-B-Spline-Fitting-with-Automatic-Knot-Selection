%% EXAMPLE: Parametric open curves
% Example 1: Heart-shaped parametric curve
% Example 2: Lissajous parametric curve

clc; 
close all; 
clear;

addpath('../../cvx')
% cvx_setup
addpath('../0. General functions and experiments')

% Selection of the example:
test_data = input('Insert the example to use (Case 1: Heart-shaped curve -- Case 2: Lissajous curve): \n');


% Reconstructed spline degree:
p = input('\nInsert the degree of the spline p (3 or 5): \n');

if p ~= 3 && p ~= 5
    error('Insert the correct degree: p=3 or p=5')
end


%% Construction of the Open Parametric Curve and sampled data:

N = 1001;  % Number of data points

switch test_data
    case 1
        % Heart-shaped parametric curve
        a=0;
        b=(3/2)*pi;

        % Function handles (for x and y coordinate)
        fx = @(t) 16*(sin(t)).^3; 
        fy = @(t) (13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t));

        % Lambda = regularization parameter
        if p == 3
            lambda = 1e-2;
        else
            lambda = 1e-6;
        end


    case 2
        % Lissajous parametric curve
        a=0;
        b=pi;

        % Function handles (for x and y coordinate)
        fx = @(t) cos(3*t);
        fy = @(t) sin(2*t);

        % Lambda = regularization parameter
        if p == 3
            lambda = 1e-2;
            
        else
            lambda = 1e-7;
        end

    otherwise
        fprintf('Invalid example number. Please enter 1 or 2.');
end

% Sampled data
t = linspace(a,b,N)';   % Parameters t
x_true = fx(t);
y_true = fy(t);
P_true = [x_true, y_true];


%% Plot original curve (Ground-truth)
figure (1)
plot(x_true, y_true, 'ro', 'MarkerFaceColor','none');
xlabel('x');
ylabel('y');
title('Original data');
grid on;
% axis equal;


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
    minimize( norm(A*C - P_true,'fro') + lambda*sum(norms(D*C,Inf,2)))
cvx_end

%% Find the active knots
jump_norm = sqrt(sum((D*C).^2,2));

idx_active = find(jump_norm > 10);    %Indices corresponding to the active knots: the knots where the jump is non zero
active_knots = inner_knots(idx_active); % Extract the active knots from the inner knots

U_active = [a*ones(1, p+1), active_knots, b*ones(1, p+1)];        % Construct the clamped knot vector with the new inner knots

fprintf('\nKnots after Step 1: %d\n', length(active_knots));
disp(active_knots);


%% Step 2:

% Define delta cluster 
delta_cluster = inner_knots(2)-inner_knots(1);

% Tolerance for the adjustment of the knots:
tol_adj = 1e-3;

% Adjust the final internal knot:
active_knots_adj = AdjustKnotsPosition(delta_cluster, tol_adj, active_knots, p, t, P_true, a, b);


fprintf('\nKnots after Step 2: %d\n', length(active_knots_adj));
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

C_final = A_final \ P_true;

% Final plot:
tt = linspace(a,b,N);
curve_final = zeros(length(tt),2);
for k = 1:length(tt)
    ik = findspan(p, tt(k), n_final, U_final_adj);
    B  = nonvanishing_basis(ik, p, tt(k), U_final_adj);
    curve_final(k,:) = B*C_final(ik-p:ik,:);
end

x_knot = fx(active_knots_adj);
y_knot = fy(active_knots_adj);
P_knot = [x_knot(:), y_knot(:)];

figure(2); 
hold on; 
grid on;
plot(P_true(:,1),P_true(:,2),'ro');
plot(curve_final(:,1),curve_final(:,2),'b-','LineWidth',2);
plot(P_knot(:,1),P_knot(:,2),'g^','LineWidth',1.5)
legend('Given data','B-spline approximation','Final knots');
title(sprintf('Final B-spline reconstruction (p = %d)',p));

%% MSE and ME
% Residual:
res = A_final * C_final - P_true;
error = vecnorm(res,2,2);  % Returns a vector of the norm 2 along the 2nd component (column),
                           % so for each row computes the 2-norm: sqrt(resx^2 + resy^2)

ME = max(error);
MSE = (norm(res,'fro')^2) / N;

fprintf('ME:  %e\n', ME);
fprintf('MSE: %e\n\n', MSE);


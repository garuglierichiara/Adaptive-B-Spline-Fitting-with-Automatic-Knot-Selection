%% Simulated examples (from Goepp 2025 sez 3.1)
clc; 
close all; 
clear;

addpath('../../cvx')
% cvx_setup

% parameters
N = 200; a = 0; b = 1;
x = linspace(a,b,N);
f = zeros (length(x),1);

% Spline degree
p = input('Insert the degree of the spline p (3 or 5): \n');

if p ~= 3 && p ~= 5
    error('Insert the correct degree: p=3 or p=5')
end

% Choose the function to test
test_case = 2; 

switch test_case
    case 1 % bump function
        for i = 1:N 
            f(i) = 0.4*(x(i)+2*exp(-(16*(x(i)-0.5))^2));
        end
  
        % lambda: parameter chosen sperimentally used 
        % in the definition of the convex optimization problem
        if p==3
            lambda = 1e-5;
        else
            lambda = 1e-8;
        end

        % PLOT TRUE FUNCTION
        figure(1);
        plot(x, f, 'r-','LineWidth',2); 
        title('BUMP FUNCTION');
        grid on;
        axis([0,1, -0.1,1.1]);

    case 2 % logit function
        for i = 1:N
            den = 1+exp(-20*(x(i)-0.5));
            f(i) = 1/den;
        end

        % lambda: parameter chosen sperimentally used 
        % in the definition of the convex optimization problem
        if p==3
            lambda = 1e-4;
        else
            lambda = 1e-7;
        end

        % PLOT TRUE FUNCTION
        figure(1);
        plot(x, f, 'r-','LineWidth',2); 
        title('LOGIT FUNCTION');
        grid on;
        axis([0,1, -0.1,1.1]);

    case 3 % sine function
        for i = 1:N
            f(i) = 0.5*sin(6*pi*x(i))+0.5;
        end

        % lambda: parameter chosen sperimentally used 
        % in the definition of the convex optimization problem
        if p==3
            lambda = 1e-5;
        else
            lambda = 1e-7;
        end

        % PLOT TRUE FUNCTION
        figure(1);
        plot(x, f, 'r-','LineWidth',2);
        title('SINE FUNCTION');
        grid on;
        axis([0,1, -0.1,1.1]);

    case 4 % spahet function
        for i = 1:N
            arg = (2*pi*(1+2^(-3/5)))/(x(i)+2^(-3/5));
            f(i) = sqrt(x(i)*(1-x(i)))*sin(arg)+0.5;
        end

        % lambda: parameter chosen sperimentally used 
        % in the definition of the convex optimization problem
        if p==3
            lambda = 1e-5;
        else
            lambda = 1e-7;
        end

        % PLOT TRUE FUNCTION
        figure(1);
        plot(x, f, 'r-','LineWidth',2);
        title('SPAHET FUNCTION');
        grid on;
        axis([0,1, -0.1,1.1]);

end

%% APPROXIMATION OF THE FUNCTION f BY A SPLINE OF DEGREE p
%% Setup parameters
             
n_init_knots = 20;                             % Initial number of knots
inner_knots = linspace(a, b, n_init_knots);    % Initial knot points
inner_knots = inner_knots(2:end-1);            % Remove the end points, consider only the inner knots

% Initial knot vector (CLAMPED)
U_init = [a*ones(1, p+1), inner_knots, b*ones(1, p+1)]; 

% RECALL: n+1 = m-p => n = m-p-1 where U = (u0=...=up, ..., um-p=...=um) (Recall in matlab you have to shift
% indices)
n_dim = length(U_init) - p - 1;
n = n_dim - 1;

%% Step 1:
%% Construction of A: A c = f
% The matrix A has dimension Nxn_dim => each row k contains the values of
% the basis fcts Ni,p evaluated in x(k)

% Pre allocation for A (sparse format)
A = sparse(N, n_dim);

% Construct A row-wise (for each element x(k))
for k = 1:N
    ik = findspan(p, x(k), n, U_init);        % Find the knot span of U_init where x(k) lives
    bv = nonvanishing_basis(ik, p, x(k), U_init);   % Find the non vanishing basis fcts (evaluated in x(k))
    % This fct returns a p+1 vector containing the p+1 non vanishing basis in
    % the knot span of x(k)
    A(k, (ik-p):ik) = bv;      % Save these values of the fcts in the k th row of A and corresponding columns (i_k-p ... i_k)
end


%% Construction of D the jump matrix
% We need to compute in each inner node the jump of the pth order
% derivative of the curve, but the pth order derivative of the curve is
% given by the lin. comb. of c times the pth derivative of the basis fcts
% Ni,p

% Consider the unique knots (= inner knots) for the jumps:
unique_inner_knots = unique(inner_knots);
num_jumps = length(unique_inner_knots);

% Pre allocate D
% Each row of D corresponds to a different internal node, while the columns
% correspond to the basis fcts
D = zeros(num_jumps, n_dim);

delta_num = 1e-4;
% Loop over each inner node for computing the jump
for j = 1:num_jumps
    u = inner_knots(j);
    % Right limit (u + delta)
    ikR = findspan(p, u + delta_num, n, U_init);
    DR = derivative_basis(ikR, p, u + delta_num, U_init, p);
    % The fct return a matrix containing in each row the derivative of
    % order k of the basis fcts non vanishing at the given knot span
    
    % Left limit (u - delta)
    ikL = findspan(p, u - delta_num, n, U_init);
    DL = derivative_basis(ikL, p, u - delta_num, U_init, p);
    % Insert the values corresponding to the derivatives of order p+1 
    % on the jth row and corresponding columns: 
    D(j, (ikR-p):ikR) = D(j, (ikR-p):ikR) + DR(p+1, :); 
    % Subtract to obtain the jump: 
    D(j, (ikL-p):ikL) = D(j, (ikL-p):ikL) - DL(p+1, :); 
end

%% Start a block to use the CVX package 
% compute c as the solution of the convex optimization problem
cvx_begin
    variable c(n_dim)      % n_dim length of the vector of the variables c
    
    % Minimize the norm of Ac-f (in this case this is a vector so we
    % use the 2 norm) + lambda * norm 1 of Dc
    minimize( norm(A*c - f, 2) + lambda * norm(D*c, 1))
cvx_end

%% 5. Find the active knots
jump = D*c;

% thereshold to select active indexes
if p == 3
    threshold = 1e-1;
else
    threshold = 1e-3;
end

idx_active = find(abs(jump)>threshold);   % Indices corresponding to the active knots: the knots where the jump is non zero
active_knots = inner_knots(idx_active);    % Extract the active knots from the inner knots

U_active = [a*ones(1, p+1), active_knots, b*ones(1, p+1)];        % Construct the clamped knot vector with the new inner knots

%% Prepare the results for the plot
% Generate 1000 points to plot the curve
xx = linspace(a, b, 1000);
curve = zeros(size(xx));

for k = 1:length(xx)
    % Use findspan and the coefficients "c" computed by cvx
    idx = findspan(p, xx(k), n, U_init);
    bv = nonvanishing_basis(idx, p, xx(k), U_init);

    % curve = linear combination of basis_fcts*coefficients c
    curve(k) = bv * c(idx-p:idx);
end

y_knots = -0.1*ones(size(active_knots)); 
% Plot
figure(2); 
hold on;
grid on;
plot(x, f, 'ro','LineWidth',1.5); 
plot(xx, curve, 'b-', 'LineWidth', 2);
% Plot the active knots as blue triangles
plot(active_knots, y_knots, 'b^', 'MarkerSize', 8,'MarkerFaceColor', 'b');
legend('Original data','Reconstructed spline', 'Active knots');
title(sprintf('Initial B-spline reconstruction and active knots (p = %d)', p));
axis([0,1, -0.3,1.1]);

%% Step 2:
% Parameter for constructing clusters chosen sperimentally
delta_cluster = inner_knots(2) - inner_knots(1);
% Tolerance for the adjustment of the knots
tol_adj = 1e-3; 
% Adjust the final internal knot:
active_knots_adj = AdjustKnotsPosition(delta_cluster, tol_adj, active_knots, p, x, f, a, b);

fprintf('\nFinal inner knots: \n')
disp(active_knots_adj);
fprintf('Number of final inner knots: %d\n', length(active_knots_adj));

%% Step 3:
% Final knot vector
U_final = [a*ones(1,p+1), active_knots_adj, b*ones(1,p+1)];
% Compute the dimension to construct A:
nF = length(U_final)-p-1;

A_final = sparse(N,nF);
for k = 1:N
    ik = findspan(p,x(k),nF-1,U_final);
    A_final(k,ik-p:ik) = nonvanishing_basis(ik,p,x(k),U_final);
end

% Solve the final Least square pb:
c_final = A_final \ f;

%% Final plot
curve_final = zeros(size(xx));

for k = 1:length(xx)
    ik = findspan(p,xx(k),nF-1,U_final);
    B  = nonvanishing_basis(ik,p,xx(k),U_final);
    curve_final(k) = B*c_final(ik-p:ik,:);
end

y_knots_adj = -0.1*ones(size(active_knots_adj)); 

figure(3); hold on; grid on;
plot(xx,curve_final,'b-','LineWidth',2);
plot(x, f, 'r', 'LineWidth',1.5); 
% Plot the final active knots as blue triangles
plot(active_knots_adj, y_knots_adj, 'b^', 'MarkerSize', 8,'MarkerFaceColor', 'b');
legend('Recovered spline','True function', 'Final knots');
title(sprintf('Final B-spline reconstruction (p = %d)',p));
axis([0,1, -0.3,1.1]);

%% MSE and ME
% Residual:
res = A_final*c_final - f;

% MSE: mean square error 
MSE = sum(res.^2) / N;

% ME: maximum error
ME = max(abs(res));

fprintf('\nMSE: %e\n', MSE);
fprintf('ME:  %e\n', ME);

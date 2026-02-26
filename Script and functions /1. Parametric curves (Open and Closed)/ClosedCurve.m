%% Example: Parametric closed (periodic) curve
% From Luo et al. (2021) - section: 5.2
clc; 
close all; 
clear;

addpath('../../cvx')
% cvx_setup
addpath('../0. General functions and experiments')

% Spline degree (for the reconstruction)
p = input('Insert the degree of the spline p (3 or 5): \n');

if p ~= 3 && p ~= 5
    error('Insert the correct degree: p=3 or p=5')
end

% End points of the domain
a = 0; 
b = 1;

% Reconstruction of the original curve (Ground truth):
p_true = 3;

% True internal knots (sec. 5.2 Luo et al. (2021))
U_true_int = [0.112, 0.325, 0.565, 0.565, 0.77, 0.832];

% Construction of a periodic knot vector (the first and last p interval
% must coincide):

% Add to the interval knots the end points of the domain:
K_true = [a, U_true_int, b]; 

% Compute the distances between the points of K:
d_true = diff(K_true);                

% Extend the knot vector to the left: sottraggo da a i passi finali (in ordine inverso)
left_ext_true = a - cumsum(d_true(end:-1:end-p_true+1));   % p values

% Extend the knot vector to the right: aggiungo a b i primi passi
right_ext_true = b + cumsum(d_true(1:p_true));             % p values

% Final knot vector (Periodic):
U_true = [left_ext_true(end:-1:1), K_true, right_ext_true];

fprintf('\nTrue inner knots: \n'); 
disp(U_true_int);

% Number of control points (total)
n_dim_tot_true = length(U_true) - p_true - 1;   % = number of control points (full)

n_tot_true = n_dim_tot_true - 1;      % n-1

% Number of dinstict control points:
% The last p control points must coincide with 
n_distinct_true = n_dim_tot_true - p_true;    % Total numer of control points (to have a spline of degree p)

fprintf('\nNumber of distinct control points:\n');
disp(n_distinct_true)

% 7 Control points DISTINCT (Obtained from the image 5.1 (b) LUO2021)
CP_dinst = [0.40    0.53;       
            0.70    1.45;
            0.90    1.30;
            0.97    0.65;
            0.50    1.85; 
            0.23    1.50;                        
            0.00   -0.50];

% Control points TOTAL (we repeat the first p points to construct a periodic Bspline of degree p) 
CP_tot = [CP_dinst; CP_dinst(1:p_true,:)];

% Sample the data points:
N = 101;                    % Number of sampled points
t = linspace(a, b, N);

% Sampled data points (from the original parametric Bspline)
P_true = zeros(N,2);
for k = 1:N
    P_true(k,:) = bspline_evaluation(p_true, t(k), n_tot_true, U_true, CP_tot); % usa la tua funzione
end

% Points to plot the curve:
tt = linspace(0,1,2001)'; 

Curve_true = zeros(length(tt),2);
for k = 1:length(tt)
    Curve_true(k,:) = bspline_evaluation(p_true, tt(k), n_tot_true, U_true, CP_tot);
end


% Plot of the original curve
figure(1); 
hold on; 
grid on; 
plot(P_true(:,1), P_true(:,2), 'ro', 'MarkerFaceColor','red');        % dati campionati (come x,f)
plot(Curve_true(:,1), Curve_true(:,2), 'r-', 'LineWidth', 2);         % curva originale
plot(CP_tot(:,1), CP_tot(:,2), 'ko-');
title('Original closed parametric curve');
legend('Sampled points P','Original closed B-spline curve', 'Control polygon');
axis([-0.02 1 -0.6 2])



%% Setup parameters:
n_init_knots = 201;     % Initial number of knots

% Construction of the initial periodic knot vector (201 points in [0,1] => 200 intervals)
n_interval = n_init_knots - 1;                 % Number of intervals
h = (b-a)/n_interval;

% Periodic extended knots (can go outside [a,b])
U_init = (a - p*h) : h : (b + p*h);

% FULL number of control points (bases) in extended representation
n_dim = n_interval + p;          
n = n_dim - 1;

% knots where we compute jumps (interior breakpoints only)
inner_knots = (a+h) : h : (b-h);    % 199 knots

% Number of inner_knots (where to compute the jump of the derivative)
num_jumps = numel(inner_knots);


%% Step 1:
%% Construction of A: A c = y_norm
A = sparse(N,n_dim);
for k = 1:N
    ik = findspan(p,t(k),n,U_init);
    A(k,ik-p:ik) = nonvanishing_basis(ik,p,t(k),U_init);
end


%% Construction of D the jump matrix
delta = 1e-2 * h;               

D = zeros(num_jumps,n_dim);

for j = 1:num_jumps
    u = inner_knots(j);

    ikR = findspan(p, u + delta, n, U_init);
    DR  = derivative_basis(ikR, p, u + delta, U_init, p);
    D(j,ikR-p:ikR) = D(j,ikR-p:ikR) + DR(p+1,:);

    ikL = findspan(p, u - delta, n, U_init);
    DL  = derivative_basis(ikL, p, u - delta, U_init, p);
    D(j,ikL-p:ikL) = D(j,ikL-p:ikL) - DL(p+1,:);
end

if p == 3
    lambda = 1e-5;          %  Regularization parameter
else % p ==5
    D = D ./ max(abs(D(:)));   % Global scaling
    % The jump operator D involves p-th order derivatives, for p=5 this makes D 
    % (and thus D*C) several orders of magnitude larger than the fitting term 
    % A*C, so the regularization can dominate and shrink the solution (curve 
    % collapse). We normalize D to make the data term and the jump-penalty 
    % comparable in scale, and then tune lambda.


    lambda = 1e-3;

end

% Start a block to use the CVX package
% In this case the variable (Control points to determine) are equal to
% n_interval, cause to impose the periodicity the last p CP must be equale
% to the first p
cvx_begin quiet
    variable Q(n_interval,2)        % Distinct control points = variables of the optimization Pb
    expression C(n_dim,2)           % Define a function of the unknown variable:
    C = [Q; Q(1:p,:)];              % copy first p points at the end (periodic)
    minimize( norm(A*C - P_true,'fro') + lambda*sum(norms(D*C,Inf,2)) )     % Objective function
cvx_end

fprintf('CP periodic check (should be close 0): %.3e\n', norm(C(1,:)-C(end-p+1,:)));



%% Find the active knots
% For a 2D Bspline curve c(t) = (cx(t), cy(t))' the jump at ti is defined
% by the maximum of the jumps from the two components x and y
% J(ti) = max{ Jx(ti), Jy(ti)}
jumps = vecnorm(D*C,Inf,2);          % num_jumps x 1
J = jumps;      % Jumps vector

eta = 5e-3;                 % Relative threshold

% Algorithm 3.1 LUO 2021
% Find the local maxima:
K0 = []; 
j = 1;
for i = 2:num_jumps-1
    jump_lmax = max([J(i-1), J(i), J(i+1)]);        % Local maximum
    if jump_lmax == J(i)     % If J(i) is a local maximum
        K0(j) = i; 
        j = j + 1;
    end
end

% Filtering using a relative threshold
K1 = find(J(K0) > eta * max(J));  % Select the indices inside K0 whose correspondent jump is above the rel threshold
idx_active = K0(K1); 

active_knots = sort(inner_knots(idx_active)); 

fprintf('\nActive knots: \n'); 
disp(active_knots);


% (Use the tt already defined before)
curve1 = zeros(length(tt),2);
for k =1:length(tt)
    ik = findspan(p, tt(k), n, U_init);
    B = nonvanishing_basis(ik, p, tt(k), U_init);
    curve1(k,:) = B*C(ik-p:ik, :);
end

% Compute the exact position (x,y) of the active knots on the reconstructed curve (After sep 1):
active_knots_final = zeros(length(active_knots), 2);
for k = 1:length(active_knots)
    % Evaluate the final spline in correspondence to the parameters U_adjusted
    active_knots_final(k,:) = bspline_evaluation(p, active_knots(k), n, U_init, C);
end

figure(3); 
hold on; 
grid on;
plot(P_true(:,1),P_true(:,2),'ro', 'MarkerFaceColor','red');
plot(Curve_true(:,1), Curve_true(:,2),'r-','LineWidth',2);
plot(curve1(:,1),curve1(:,2),'b-','LineWidth',1.5);
plot(active_knots_final(:,1), active_knots_final(:,2), 'k^', 'MarkerSize', 8, 'Linewidth',1.8,'MarkerFaceColor', 'g');
legend('Original data points', 'Original B-spline','Reconstructed spline', 'Active knots');
title(sprintf('Initial B-spline reconstruction and active knots (p = %d)', p));



%% Step 2:

% Define delta cluster:
h = inner_knots(2) - inner_knots(1);

% Tune the parameters depending on p:
if p == 3
    delta_cluster = h;        
    % Tolerance for the adjustment of the knots:
    tol_adj = 1e-5;
else
    delta_cluster = 3*h;
    tol_adj = 1e-7;
end


% Adjust the final internal knot:
active_knots_adj = AdjustKnotsPosition_closed(delta_cluster, tol_adj, active_knots, p, t, P_true, a, b);

% Remove the knots outside (a,b) => internal knots
U_adjusted = active_knots_adj(active_knots_adj>a & active_knots_adj<b);

fprintf('\nFinal inner knots: \n'); 
disp(U_adjusted);
fprintf('Number of final inner knots: %d\n', length(U_adjusted));

%% Final spline reconstruction:
K = [a, U_adjusted, b];
dK = diff(K);

% Estensione periodica a sinistra/destra basata sugli intervalli di K
left_ext  = a - cumsum(dK(end:-1:end-p+1));
left_ext  = left_ext(end:-1:1);
right_ext = b + cumsum(dK(1:p));

U_final = [left_ext, K, right_ext];

n_dim_final = length(U_final) - p - 1;   % final number of  Control points (total)
n_final = n_dim_final - 1;

A_final_full = sparse(N, n_dim_final);
for k = 1:N
    ik = findspan(p, t(k), n_final, U_final);
    A_final_full(k, ik-p:ik) = nonvanishing_basis(ik, p, t(k), U_final);
end



% In a periodic B_spline, some control points are repeated, they are not
% new: indeed the last pa are copied from the firs p points.
% Thus, we cannot consider all the control points as unknown independent 
% variables in the LS Pb, otherwise we would get a double estimate of the
% same points. 
% So we solve the LS only for the distinct points Q and the we constuct
% C = [Q; Q(1:p)]

% How can we obtain AF?
% Start from the eq. A_full C_full = P_true
% Since C_full = [Q; Q(1:p)] = > A_full C_full = A1 Q + A2 Q(:,1:p)
% with A1 = A_full(:,1:M2) e A2 = A_full(:, M2+1:M2+p)
% So the coefficients of Q(:,1:p) are given by A1 + A2 while the coeff. of the
% other columns only by A1

M2 = n_dim_final - p;       % Number of distinct CP
AF = A_final_full(:, 1:M2); % Initialize AF = A1 (fist M2 columns of A_full)

% Loop over the first p columns of AF (for these we need to sum A2):
for r = 1:p
    AF(:, r) = AF(:, r) + A_final_full(:, M2 + r);
end

% LS for the distinct control points:
QF = AF \ P_true;            % (M2 x 2)
CF = [QF; QF(1:p,:)];        % full periodic CP


% (Use the previous tt)
curve_final = zeros(length(tt),2);

for k = 1:length(tt)
    ik = findspan(p,tt(k),n_final,U_final);
    B  = nonvanishing_basis(ik,p,tt(k),U_final);
    curve_final(k,:) = B*CF(ik-p:ik,:);
end

% Compute the exact position (x,y) of the knots on the reconstructed curve:
knots_final = zeros(length(U_adjusted), 2);
for k = 1:length(U_adjusted)
    % Evaluate the final spline in correspondence to the parameters U_adjusted
    knots_final(k,:) = bspline_evaluation(p, U_adjusted(k), n_final, U_final, CF);
end

figure(4); 
hold on; 
grid on;
plot(P_true(:,1),P_true(:,2),'ro', 'MarkerFaceColor','red');
plot(curve_final(:,1),curve_final(:,2),'b-','LineWidth',2.5);
plot(Curve_true(:,1),Curve_true(:,2),'r-','LineWidth',1);
plot(knots_final(:,1), knots_final(:,2), 'k^', 'MarkerSize', 8, 'Linewidth',1.8,'MarkerFaceColor', 'g');
legend('Original data points', 'Reconstructed spline', 'Original B-spline', 'Final knots');
title(sprintf('Final B-spline reconstruction (p = %d)',p));


%% MSE and ME
% Residual:
res = A_final_full*CF - P_true;
error = vecnorm(res,2,2);  % Returns a vector of the norm 2 along the 2nd component (column),
                           % so for each row computes the 2-norm: sqrt(resx^2 + resy^2)

% MSE: mean square error 
MSE = (norm(res,'fro')^2) / N;

% ME: maximum error
ME = max(error);

fprintf('\nMSE: %e\n', MSE);
fprintf('ME:  %e\n', ME);
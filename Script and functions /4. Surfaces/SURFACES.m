%% Example s_spline.db taken from LAB 7 material
% IN THIS SCRIPT WE EXTEND AND TEST THE KANG'S ALGORITHM (STEP 1: SPARSE OPTIMIZATION
% AND ACTIVE KNOTS SELECTION + STEP 2: ADJUSTMENT AND PRUNING) ON A
% SURFACE.
% -----------------------------------------------------------------------

% Upload files
clc; clear; close all;

addpath('../../cvx')
% cvx_setup
addpath('../0. General functions and experiments')

%% -------------------------------------------------------------
filename = 's_spline.db'; 
surfspline = surf_spline_load(filename);

%-------------------------------------------------
% 1. Extract the parameters from the structure.
% ------------------------------------------------
p_true = surfspline.deguv(1);
q_true = surfspline.deguv(2);
P = surfspline.cp;
U = surfspline.ku;
V = surfspline.kv;

[nu, mv, ~] = size(P);
n = nu - 1; 
m = mv - 1;

% -------------------------------------------------
% 2.  Surface
% -------------------------------------------------
N = 50;     % Number of sampled
u_vec = linspace(U(p_true+1), U(end-p_true), N);  % evaluation points in the u direction
v_vec = linspace(V(q_true+1), V(end-q_true), N);  % evaluation points in the v direction

% initialization of the surface components
Sx = zeros(N);
Sy = zeros(N);
Sz = zeros(N);

% Evaluate the surface for plotting it.
for k = 1:N
    for l = 1:N
        S = Bspline_surface_eval([u_vec(k) v_vec(l)], P, p_true, q_true, U, V);
        Sx(k,l) = S(1); 
        Sy(k,l) = S(2); 
        Sz(k,l) = S(3);
    end
end

% ------------------------------------------------------------
% 3. Plot the true surface and the control net.
% ------------------------------------------------------------
figure(1); 
hold on; 
grid on;
% Surface
surf(Sx, Sy, Sz, 'FaceAlpha', 0.8, 'EdgeColor', 'none'); 
 
% Plot Control Net 
for i = 1:nu
    plot3(P(i,:,1), P(i,:,2), P(i,:,3), 'k-o');
end
for j = 1:mv
    plot3(P(:,j,1), P(:,j,2), P(:,j,3), 'k-o');
end

view(3); 
title(sprintf('B-spline Surface (p=%d, q=%d)', p_true, q_true));


%% ==============================================================
% 4. KANG APPLICATION
% ===============================================================

% Spline degree:
p = 3;
q = 3;

% Parameters
n_init_u = 40;       % Number initial knots in u 
n_init_v = 40;       % Number initial knots in v
lambda = 0.1;        

% Exstreme points definition
ua = U(p_true+1);         
ub = U(end-p_true);
va = V(q_true+1); 
vb = V(end-q_true);

% Initial inner knots definition
inner_u = linspace(ua, ub, n_init_u); 
inner_u = inner_u(2:end-1);
inner_v = linspace(va, vb, n_init_v); 
inner_v = inner_v(2:end-1);


% Initial Knot Vectors (CLAMPED)
U_init = [ua*ones(1,p+1), inner_u, ub*ones(1,p+1)];
V_init = [va*ones(1,q+1), inner_v, vb*ones(1,q+1)];


n_dimU = length(U_init) - p - 1;
n_dimV = length(V_init) - q - 1;
n_u_opt = n_dimU - 1;
n_v_opt = n_dimV - 1;
total_vars = n_dimU * n_dimV;

% Construcion matrix A (Fitting)
% We transform the surface data 2D in a fitting problem in 1D (vector)
X_target = Sx(:); 
Y_target = Sy(:); 
Z_target = Sz(:);
Target = [X_target, Y_target, Z_target];
N_tot = size(Target,1);

A = sparse(N_tot, total_vars);  % it is sparse because for each (u,v) just 
% (p+1)(q+1) basis are not zero → so for each row we have few non zero
% entries

% Loop over the points (u_vec, v_vec) 
for k = 1:N %u
    for l = 1:N  %v
        row = (l-1)* N + k;       % To save in the correct position the values of A = Ni(uk)Mj(vl)
        % The '3D vector' Target contains uv1,...,uvn
        % => for the point (u(k),v(l)) we have to skip the previous rows
        % ((l-1) blocks of N elements and then in the l-th block we
        % consider the kth row => (l-1)*N + k
        u = u_vec(k); 
        v = v_vec(l);
        
        iu = findspan(p, u, n_u_opt, U_init); 
        iv = findspan(q, v, n_v_opt, V_init); 
        
        Nu = nonvanishing_basis(iu, p, u, U_init);
        Nv = nonvanishing_basis(iv, q, v, V_init);
        
        % Loop over the control points
        for jj = 0:q
            for ii = 0:p
                col_idx = (iv-q+jj-1)*n_dimU + (iu-p+ii);
                A(row, col_idx) = Nu(ii+1) * Nv(jj+1);
            end
        end
    end
end

%% STEP 1: SPARSE FITTING
% 1. Construction of the jump's matrices Du and Dv
delta = 1e-6;

% Jump in U
Du_base = zeros(length(inner_u), n_dimU);
for j = 1:length(inner_u)
    u_node = inner_u(j);
    % Right limit (u + delta)
    ikR = findspan(p, u_node + delta, n_u_opt, U_init); 
    DNR = derivative_basis(ikR, p, u_node + delta, U_init, p);
    % The fct return a matrix containing in each row the derivative of
    % order k of the basis fcts non vanishing at the given knot span
    
    % Insert the values corresponding to the derivatives of order p+1 
    % on the jth row and corresponding columns: 
    Du_base(j, (ikR-p):ikR) = DNR(p+1, :);
    % Left limit (u - delta)
    ikL = findspan(p, u_node - delta, n_u_opt, U_init); 
    DNL = derivative_basis(ikL, p, u_node - delta, U_init, p);
    % Subtract to obtain the jump: 
    Du_base(j, (ikL-p):ikL) = Du_base(j, (ikL-p):ikL) - DNL(p+1, :);
end
% Traasform the jump matrix in 1D into 2D
Du = kron(eye(n_dimV), Du_base);   % Jump deriv of order p

% Jump in V
Dv_base = zeros(length(inner_v), n_dimV);
for j = 1:length(inner_v)
    v_node = inner_v(j);
    % Right limit (v + delta)
    ikR = findspan(q, v_node + delta, n_v_opt, V_init); 
    DNR = derivative_basis(ikR, q, v_node + delta, V_init, q);
    % The fct return a matrix containing in each row the derivative of
    % order k of the basis fcts non vanishing at the given knot span
    
    % Insert the values corresponding to the derivatives of order p+1 
    % on the jth row and corresponding columns: 
    Dv_base(j, (ikR-q):ikR) = DNR(q+1, :);
    % Left limit (v - delta)
    ikL = findspan(q, v_node - delta, n_v_opt, V_init); 
    DNL = derivative_basis(ikL, q, v_node - delta, V_init, q);
    % Subtract to obtain the jump: 
    Dv_base(j, (ikL-q):ikL) = Dv_base(j, (ikL-q):ikL) - DNL(q+1, :);
end
Dv = kron(Dv_base, eye(n_dimU));

% Normalization
Du = Du ./ max(abs(Du(:)));
Dv = Dv ./ max(abs(Dv(:)));

% 2. Optimization CVX
cvx_begin
    variable C(total_vars, 3) 
    minimize( norm(A*C - Target, 'fro') + lambda * ( sum(norms(Du*C,Inf,2)) + sum(norms(Dv*C,Inf,2))) )
    % minimize( norm(A*C - Target, 'fro') + lambda * (norm(Du*C, 1) + norm(Dv*C, 1)) )        
cvx_end

% 3. Active knots selection
jump_u = Du * C;  % It contains the jumps of the p-th derivative in u across each u=uj internal node and for each v-fixed strip.      
jump_v = Dv * C;  % It contains the jumps of the q-th derivative in v across each v=vj internal node and for each u-fixed strip.

% for each node, compute the jump norm
norm_ju = vecnorm(jump_u,Inf,2);
norm_jv = vecnorm(jump_v,Inf,2);

% Select active indexes for inner_u e inner_v
mat_ju = reshape(norm_ju, length(inner_u), n_dimV);
% mat_ju(i,j) = jump magnitude at u = inner_u(i), for the j-th v-basis strip
mean_ju = mean(mat_ju, 2);   % average jump over all v-strips for each inner_u(i)

mat_jv = reshape(norm_jv, n_dimU, length(inner_v));
mean_jv = mean(mat_jv, 1);

idx_u_active = find(mean_ju > 1e-4);  
idx_v_active = find(mean_jv > 2e-4);  

inner_u_final = inner_u(idx_u_active);
inner_v_final = inner_v(idx_v_active);

U_final = [ua*ones(1,p+1), inner_u_final, ub*ones(1,p+1)];
V_final = [va*ones(1,q+1), inner_v_final, vb*ones(1,q+1)];

fprintf('\nStep 1 completed.\nActive knots U: %d\nActive knots V: %d\n', ...
    length(inner_u_final), length(inner_v_final));


figure(2); 
hold on; 
grid on; 
axis equal;
[UU_k, VV_k] = meshgrid(inner_u_final, inner_v_final);
plot(UU_k(:), VV_k(:), 'b^', 'MarkerFaceColor','b');
xlabel('u'); ylabel('v');
title('Active knots grid');



%% STEP 2: ADJUSTMENT
% delta_cluster for inner_u_final and for inner_v_final, and tollerance
delta_u = 2e-2; 
delta_v = 0.05; 
tol_adj = 1e-3;

% 1. Adjust active knots in U using active inner knots in v
inner_u_adj = AdjustKnotsPosition_Surf(delta_u, tol_adj, inner_u_final, p, inner_v_final, q, u_vec, v_vec, Sx, Sy, Sz, 'u', ua, ub, va, vb);

% 2. Adjust active knots in V using the adjusted knots in U
inner_v_adj = AdjustKnotsPosition_Surf(delta_v, tol_adj, inner_v_final, q, inner_u_adj, p, u_vec, v_vec, Sx, Sy, Sz, 'v', ua, ub, va, vb);

U_final_adj = [ua*ones(1,p+1), inner_u_adj, ub*ones(1,p+1)];
V_final_adj = [va*ones(1,q+1), inner_v_adj, vb*ones(1,q+1)];

fprintf('\nStep 2 completed.\nFinal knots U: %d\nFinal knots V: %d\n', ...
    length(inner_u_adj), length(inner_v_adj));

%% STEP 3: FINAL FITTING AND RECONSTRUCTION
n_dimU_f = length(U_final_adj) - p - 1;
n_dimV_f = length(V_final_adj) - q - 1;
total_vars_f = n_dimU_f * n_dimV_f;

% 1. A_final matrix construction
A_final = sparse(N_tot, total_vars_f);
for k = 1:N
    for l = 1:N
        row = (l-1)*N + k;
        iu = findspan(p, u_vec(k), n_dimU_f-1, U_final_adj);
        iv = findspan(q, v_vec(l), n_dimV_f-1, V_final_adj);
        Nu = nonvanishing_basis(iu, p, u_vec(k), U_final_adj);
        Nv = nonvanishing_basis(iv, q, v_vec(l), V_final_adj);
        for jj = 0:q
            for ii = 0:p
                col_idx = (iv-q+jj-1)*n_dimU_f + (iu-p+ii);
                A_final(row, col_idx) = Nu(ii+1) * Nv(jj+1);
            end
        end
    end
end

% 2. Resolution of the least square problem for the control points C
% We solve separately the three coordinates X, Y, Z
C_final = A_final \ [Sx(:), Sy(:), Sz(:)]; 

% 3. Valutation reconstructed surface for the plot
[UU, VV] = meshgrid(linspace(ua,ub,50), linspace(va,vb,50));
Sx_rec = zeros(50); 
Sy_rec = zeros(50); 
Sz_rec = zeros(50);

% Transform C_final into an array 3D for the valutation
P_final = reshape(C_final, [n_dimU_f, n_dimV_f, 3]);

for k = 1:50
    for l = 1:50
        S = Bspline_surface_eval([UU(k,l) VV(k,l)], P_final, p, q, U_final_adj, V_final_adj);
        Sx_rec(k,l) = S(1); 
        Sy_rec(k,l) = S(2); 
        Sz_rec(k,l) = S(3);
    end
end

figure(3); 
s1 = surf(Sx_rec, Sy_rec, Sz_rec, 'FaceAlpha', 0.9, 'EdgeColor', 'none'); 
hold on; 
colormap jet; 
camlight; 
lighting phong;
grid on;
view(3);       
axis tight;   
title(sprintf('Reconstructed surface (p=%d, q=%d)', p, q));

figure(4); 
hold on; 
grid on; 
axis equal;
[UU_k, VV_k] = meshgrid(inner_u_adj, inner_v_adj);
plot(UU_k(:), VV_k(:), 'b^', 'MarkerFaceColor','b');
xlabel('u'); ylabel('v');
title('Final knot grid');



% Residual
res = A_final * C_final - Target;

% Pointwise Euclidean error
err = vecnorm(res, 2, 2);   % sqrt(dx^2 + dy^2 + dz^2)

% MSE: mean square error 
MSE = norm(res,'fro')^2/N;

% ME: maximum error
ME  = max(err);

fprintf('\nMSE: %e\n', MSE);
fprintf('ME:  %e\n', ME);

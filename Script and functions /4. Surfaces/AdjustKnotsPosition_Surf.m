%% THE FUNCTION AdjustKnotsPosition_Surf ADJUSTS THE POSITION OF THE ACTIVE KNOTS IN U_TILDE WHILE
% FREEZING THE ACTIVE KNOTS IN OTHER_KNOTS

function [UU_adj] = AdjustKnotsPosition_Surf(delta, tol, U_tilde, p, other_knots, other_p, u_vec, v_vec, Sx, Sy, Sz, direction, ua, ub, va, vb)
% INPUTS:
% - delta: knot spacing
% - tol: threshold for stopping the bisection method
% - U_tilde: active knots (inner points) to be adjusted
% - p: surface degree in the 'direction'
% - other_knots: active knots in the other direction to keep freezed
% - other_p: surface degree in the other direction
% - u_vec: evaluation points in the u direction
% - v_vec: evaluation points in the v direction
% - Sx, Sy, Sz: surface components
% - direction: 'u' if we are optimizing inner_u, 'v' if we are optimizing inner_v
% - ua, ub, va, vb: domain end-points
% ---------------------------------------------
% OUTPUT:
% - UU_adj: final knots adjusted in the direction 'direction' (inner knots)
% ---------------------------------------------

% Cluster classification:
G = ClusterClassify(delta, U_tilde);
UU = sort(U_tilde(G)); % Extract and sort the initial and end knots of each cluster
% UU = [s1 e1 s2 e2 ...] where si: knot where the i-th cluster starts
%                              ei: knot where the i-th cluster ends
% Force UU to be a row vector
UU = UU(:)'; 

i = 1;        % initialization of index i

% Loop over each cluster [UU(i), UU(i+1)]
while (i < length(UU))
    pre_val = UU(i);
    next_val = UU(i+1);
    % 1. Decide if the interval should be represented by a single knot, or 
    % a double knot at the same location.
    % The decision is based on comparing LS errors.

    % Compute the mid point:
    mid = 0.5 * (pre_val + next_val);
    
    % Mdt: inserting a single knot at the midpoint (Merge)
    Mdt = [UU(1:i), mid, UU(i+2:end)];
    % Mpt: inserting a double knot at the midpoint (Split, increase the multiplicity)
    Mpt = [UU(1:i), mid, mid, UU(i+2:end)];

    % Compute the LS error in the 'direction'
    if strcmp(direction, 'u')
        % the first argument are knots in the u direction and it is the one
        % that changes. The second argument is the one in the v direction
        % and it is fixed.
        Ed = LeastSquareFit_Surf(Mdt, other_knots, p, other_p, u_vec, v_vec, Sx, Sy, Sz, ua, ub, va, vb);
        Ep = LeastSquareFit_Surf(Mpt, other_knots, p, other_p, u_vec, v_vec, Sx, Sy, Sz, ua, ub, va, vb);
    else
        % the second argument are knots in the v direction and it is the one
        % that changes. The first argument is the one in the u direction
        % and it is fixed.
        Ed = LeastSquareFit_Surf(other_knots, Mdt, other_p, p, u_vec, v_vec, Sx, Sy, Sz, ua, ub, va, vb);
        Ep = LeastSquareFit_Surf(other_knots, Mpt, other_p, p, u_vec, v_vec, Sx, Sy, Sz, ua, ub, va, vb);
    end

    % Compare the errors
    if 2 * Ep < Ed
        flag = 1; % Split case (double knot) 
    else
        flag = 2; % Merge case (single knot)
    end
    
    % 2. Adjustment strategy (Bisection method)
    
    % Initialize the starting points
    temp_pre = pre_val;
    temp_next = next_val;

    while (temp_next - temp_pre > tol)    % Stopping criterion: we stop the bisection when the interval is small enough
        mid_point = 0.5 * (temp_next + temp_pre);
        
        % Try to move the left boundary knot  
        Lt = [UU(1:i), mid_point, UU(i+2:end)];
        % Try to move the right boundary knot  
        Rt = [UU(1:i-1), mid_point, UU(i+1:end)];
        
        % Compute the LS error in the 'direction'
        if strcmp(direction, 'u')
            % the first argument are knots in the u direction and it is the one
            % that changes. The second argument is the one in the v direction
            % and it is fixed.
            Lerr = LeastSquareFit_Surf(Lt, other_knots, p, other_p, u_vec, v_vec, Sx, Sy, Sz, ua, ub, va, vb);
            Rerr = LeastSquareFit_Surf(Rt, other_knots, p, other_p, u_vec, v_vec, Sx, Sy, Sz, ua, ub, va, vb);
        else
            % the second argument are knots in the v direction and it is the one
            % that changes. The first argument is the one in the u direction
            % and it is fixed.
            Lerr = LeastSquareFit_Surf(other_knots, Lt, other_p, p, u_vec, v_vec, Sx, Sy, Sz, ua, ub, va, vb);
            Rerr = LeastSquareFit_Surf(other_knots, Rt, other_p, p, u_vec, v_vec, Sx, Sy, Sz, ua, ub, va, vb);
        end
        
        % Compare the errors and move the boundary returning a low error
        if Rerr < Lerr
            temp_pre = mid_point;    % move search window to the right half
        else
            temp_next = mid_point;   % move search window to the left half
        end
    end
    
    % Update the final vector UU
    final_pos = 0.5 * (temp_pre + temp_next);  % Compute the final mid point (Adjusted knot)
    
    % Choose the multiplicity (based on the previous decision)
    if flag == 2 
        % Merge: The interval is substituted by a unique knot 
        UU(i) = final_pos;
        UU(i+1) = [];
        % Do not increment i here because the length of the vector is
        % reduced
    else
        % Double: We substitute the end points of the interval with the
        % "new" knot with multiplicity = 2
        UU(i) = final_pos;
        UU(i+1) = final_pos;
        i = i + 1;      % Increment i to skip the second repeated knot
    end
    i = i + 1;          % Pass to the following "cluster"
end
UU_adj = UU;
end

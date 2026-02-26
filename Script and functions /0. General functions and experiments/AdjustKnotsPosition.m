%% THE FUNCTION AdjustKnotsPosition ADJUSTS THE POSITION OF THE ACTIVE KNOTS
%
% It follows the algorithm 4: Locally adjust knots position (II) from 
% Kang et al. (2015)
% The function works with non-parametric and parametric curves (open).

function [UU] = AdjustKnotsPosition(delta, tol, U_tilde, p, x, y, a, b)
% INPUTS:
% - delta: knot spacing
% - tol: threshold for stopping the bisection method
% - U_tilde: active knots (inner points)
% - p: spline degree
% - x: parameters (Nx1)
% - y: data values -> (Nx1) for scalar functions or (NxD) for parametric curves
% - a,b: domain end-points
% ---------------------------------------------
% OUTPUT:
% - UU: final knots adjusted (inner knots)
% ---------------------------------------------

% Force x to be a column vector
x = x(:);

% If y is a vector => force it to be a column vec.
if (size(y,1)==1 || size(y,2) == 1)      
    y = y(:);                
end

% Cluster classification:
G = ClusterClassify(delta, U_tilde);
UU = U_tilde(G);        % Extract the initial and end knots of each cluster
% UU = [s1 e1 s2 e2 ...] where si: knot where the i-th cluster starts
%                              ei: knot where the i-th cluster ends
UU = sort(UU(:))';         % Force UU to be a row vector ordered


i = 1;  % initialization of index i

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
    Mdt = [UU(1:i), mid, UU(i+1:end)];
    Ed = LeastSquareFit(Mdt, p, x, y, a, b);    % Compute the LS error

    % Mpt: inserting a double knot at the midpoint (Split, increase the multiplicity)
    Mpt = [UU(1:i), mid, mid, UU(i+1:end)];
    Ep = LeastSquareFit(Mpt, p, x, y, a ,b);    % Compute the LS error

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

    while (temp_next - temp_pre > tol)      % Stopping criterion: we stop the bisection when the interval is small enough
        
        mid_point = 0.5 * (temp_next + temp_pre);   % Mid point

        % Try to move the left boundary knot  
        Lt = [UU(1:i), mid_point, UU(i+2:end)];
        Lerr = LeastSquareFit(Lt, p, x, y, a ,b);   % Compute the error

        % Try to move the right boundary knot  
        Rt = [UU(1:i-1), mid_point, UU(i+1:end)];
        Rerr = LeastSquareFit(Rt, p, x, y, a ,b);   % Compute the error

        % Compare the errors and move the boundary returning a low error
        if Rerr < Lerr
            temp_pre = mid_point; % move search window to the right half
        else
            temp_next = mid_point; % move search window to the left half
        end
    end

    % 3. Update the final vector UU
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
        i = i+1;    % Increment i to skip the second repeated knot
    end
    i = i + 1;      % Pass to the following "cluster"
end
end

%% THE FUNCTION fittness_function constructs the objective function used by Differential Evolution (DE)
% Given a candidate set of INNER KNOTS, this function:
%   1) builds the full clamped knot vector
%   2) builds the collocation matrix A
%   3) solves a least-squares problem to estimate control points
%   4) returns the mean squared fitting error (MSE).
%
% DE calls this function many times. The goal of DE is to MINIMIZE the MSE
% by moving the knot locations.

function mse = fitness_function(settings, params)
% INPUTS:
% - settings: structure that contains all the fixed data that DE passes to
% fittness (as the data points, spline degree...)
% - params: parameters optimized by DE (inner knot points)
% ---------------------------------------------
% OUTPUT:
% - mse: mean square error (scalar)
% ---------------------------------------------

% DE can pass the parametes as a structure or as a vector
% We transform them alwasy in a vector (if they aren't)
% Trasforma params in vettore se è una struct
    
    if isstruct(params)
        p = cell2mat(struct2cell(params)); 
    else
        p = params; 
    end
    
    % Save the values from settings:
    d = settings.deg;       % Spline degree
    P = settings.P;         % Data Points
    t = settings.t;         % Parameter values (associated to the data points)
    a = settings.a;         % Domain end points
    b = settings.b;
    
    % Candidate inner knots (sorted)
    U_inner = sort(p(:))';
    
    % Erase invalid candidates (knots too close):
    % If two knots are almost equal, the spline basis becomes ill-conditioned
    % and the least-squares matrix A can be nearly singular.
    % In that case we return a very large error to penalize the candidate.
    if any(diff(U_inner) < 1e-8)
        mse = 1e12; 
        return; 
    end

    % Construct the full clamped vector:
    U_full = [a*ones(1, d+1), U_inner, b*ones(1, d+1)];
    n_bases = length(U_full) - d - 1;
    N = size(P, 1);

    % We use "try and catch" because we could have problems in computing
    % the MSE => so if everything works well we get from "try" the mse
    % computed, otherwise we fix in "catch" a high value => so DE
    % understands that this solution is not good at all
    try
        % Build the collocation matrix:
        A = zeros(N, n_bases);
        for k = 1:N
            ik = findspan(d, t(k), n_bases-1, U_full);
            A(k, ik-d:ik) = nonvanishing_basis(ik, d, t(k), U_full);
        end
        
        % We use the pseudo inverse  pinv(A) instead of A\P to be more 
        % robust when A is close to singular 
        C_temp = pinv(A) * P; 
        
        % Compute the mse
        res = A * C_temp - P;
        mse = mean(sum(res.^2, 2));
    catch
        mse = 1e10;
    end
end

%% THE FUNCTION bspline_evaluation EVALUATES THE B-SPLINE CURVE C IN u

function [C_u] = bspline_evaluation(p, u, n, U, P)
% INPUTS:
% - p: degree
% - u: evaluation point
% - n: n+1 is the dimension of spline space
% - U: clamped knot vector
% - P: Matrix of control points
% -----------------------------------------
% OUTPUT:
% - C_u: vector corresponding to the value of the B-spline curve C evaluated in u
% -----------------------------------------

i = findspan(p, u, n, U);
N = nonvanishing_basis(i, p, u, U);

% Force P to have the control points as rows:
dim = min(size(P));
if size(P,2)~= dim
    P = P';
end

% Pre allocate the variacle
C_u = zeros(1, dim);

% Loop over the non vanishing basis at u
for j = 0 : p
    basis_val = N(j+1);             % N contains only the p+1 basis elements non zero
    control_point = P(i-p+j, :);    % Corresponding control point
    C_u = C_u + basis_val * control_point;  % Update the sum
end
end
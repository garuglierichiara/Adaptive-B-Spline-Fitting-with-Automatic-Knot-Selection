%% THE FUNCTION Bspline_surface_eval EVALUATE THE B-SPLINE SURFACE S IN A POINT x
%
% This function evaluates a B-spline surface
%   S(u,v) = sum_{i=0}^{n} sum_{j=0}^{m} P_{i,j} N_{i,p}(u) N_{j,q}(v)
% at a given parameter pair (u,v).

function [S] = Bspline_surface_eval(x, P, p, q, U, V)
% INPUTS:
% - x: evaluation point, vector of 2 components x = [u,v]
% - P: control points net, size (n+1) x (m+1) x d
% - p: degree on u direction
% - q: degree on v direction
% - U: clamped knot vector in u direction
% - V: clamped knot vector in v direction
% -----------------------------------------
% OUTPUT:
% - S: vector corresponding to the value of the B-spline surface S
% evaluated in x, size (1xd)
% -----------------------------------------

% Extract the u and v coordinates:
u = x(1); 
v = x(2);

% Last index in control net (u and v direction respectively):
n = size(P,1) - 1;
m = size(P,2) - 1;

% Find the knot spans containing u and v:
i = findspan(p, u, n, U);
j = findspan(q, v, m, V);

% Compute the non vanishing basis functions at (u,v)
Nu = nonvanishing_basis(i, p, u, U);
Nv = nonvanishing_basis(j, q, v, V);

% Pre allocate the vector that will contain the value S(u,v)
d = size(P,3);
S = zeros(1,d);

% Loop over v active basis functions (q+1 terms)
for l = 0:q
    temp = zeros(1,d);  % Initialize the inner sum over u  = 0

    for k = 0:p     % Loop over u active basis functions (p+1 terms)
        % Compute the correct indices for the control net
        cp_u = i - p + k;
        cp_v = j - q + l;
        
        % Extract the control point coordinates P(cp_u, cp_v, :) and
        % reshape in a vector (1xd)
        cp_coord = reshape(P(cp_u, cp_v, :), 1, d);
        
        temp = temp + Nu(k+1) * cp_coord;       % Update the inner sum over u
    end
    % Multiply by the v basis and add the term to the final sum (S)
    S = S + Nv(l+1) * temp;
end
end
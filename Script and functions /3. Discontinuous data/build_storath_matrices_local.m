%% THE FUNCTION build_storath_matrices_local CONSTRUCTS THE MATRIX Ar AND THE VECTOR yr 
% USED TO SOLVE THE FIRST OPTIMIZATION PROBLEM OF THE STORATH AND WEINMANN STRATEGY TO
% FOUND THE SET OF DISCONTINUITIES OF A DISCONTINOUS FUNCTION f.
% ------------------------------------------------

function [Ar, yr] = build_storath_matrices_local(x, y, p, noise_level)
% INPUTS:
% - x: vector of evaluation points
% - y: vector of noisy data
% - p: stiffness parameter
% - noise_level: standard deviation
% -------------------------------------------------
% OUTPUTS:
% - Ar: (3r-2) x (2r) matrix 
% - yr: (3r-2) column vector
% where r = length(x)
% -------------------------------------------------

r = length(x);

% Parameters defined to construct the matrix Ar and the vector yr
alfa = sqrt(p) / noise_level;
beta = sqrt(1 - p);

% Initialization
nrows = 3*r - 2;
Ar = zeros(nrows, 2*r);
yr = zeros(nrows, 1);

% Fidelity blocks
for i = 1:r
    Ar(3*i-2, 2*i-1) = alfa;
    yr(3*i-2) = alfa * y(i);
end

% Regularization blocks
for i = 1:r-1
    di = x(i+1) - x(i);
    Vi = [ 2*sqrt(3)/di^(1.5),  sqrt(3)/sqrt(di); 0, 1/sqrt(di) ];
    Wi = [ -2*sqrt(3)/di^(1.5), sqrt(3)/sqrt(di); 0, -1/sqrt(di) ];

    Ar(3*i-1:3*i, 2*i-1:2*i) = beta * Vi;
    Ar(3*i-1:3*i, 2*i+1:2*i+2) = beta * Wi;
end
end
%% THE FUNCTION nonvanishing_basis COMPUTES THE NON-VANISHING BASIS FUNCTIONS OF DEGREE p EVALUATED
% AT u BELONGING TO [U(i), U(i+1))

function N = nonvanishing_basis(i,p,u,U)
% INPUTS:
% - i: knot span index which contains u
% - p: degree
% - u: evaluation points
% - U: clamped knot vector
% -----------------------------------------
% OUTPUT:
% - N: array N = (Ni-p,p(u),...,Ni,p(u))
% where its compontens are the non-vanishing basis functions of degree p
% evaluated at the point u.
% -----------------------------------------

N = zeros(1,p+1);
left = zeros(1,p);
right = zeros(1,p);
N(1) = 1;           % Ni,0(u) = 1         
        
for j=1:p   % Loop over the degree
    left(j) = u - U(i+1-j);         % Num. for column 1
    right(j)= U(i+j) - u;           % NUm. for column 2
    saved = 0;
    for r = 1:j   % Loop to compute the non -vanishing basis of degree j (j+1 elements)   
        denom = right(r) + left(j-r+1);
        if denom == 0
            temp = 0;
        else
            temp = N(r) / denom;
        end
        N(r) = saved + right(r)*temp;
        saved = left(j-r+1)*temp;
    end
    N(j+1)=saved;
end
end
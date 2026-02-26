%% THE FUNCTION findspan DETERMINES i s.t. u BELONGS TO [U(i), U(i+1))

function i = findspan(p, u, n, U)
% INPUTS:
% - p: degree
% - u: evaluation point
% - n: dimension of the spline space
% - U: knot vector
% ---------------------------------------------
% OUTPUT:
% - i: index in which the evaluation point u belongs
% ---------------------------------------------

if u >= U(n+2)
    i = n+1; 
    return;
end
if u <= U(p+1)
    i = p+1; 
    return;
end

low = p+1;
high = n+2;
mid = floor((low+high)/2);

while high-low > 1
    if u < U(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low+high)/2);
end
i = low;
end
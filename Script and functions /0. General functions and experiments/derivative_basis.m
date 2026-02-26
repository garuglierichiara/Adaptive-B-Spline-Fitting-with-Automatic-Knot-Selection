%% THE FUNCTION derivative_basis EVALUATE THE DERIVATIVES OF THE NONZERO
%% B-SPLINE BASIS AT u IN [U(i), U(i+1))

% Formula:
% N'_i,p(u) = p* [1/(u_i+p - ui) * N_i,p-1 + 1/(u_i+1 - u_i+p+1) * N_i+1,p-1]

% => N^(k)_i,p(u) = p!/(p-k)  sum_{i=0,..,k} ( a_k,j * N_i+j,p-k )

function DN = derivative_basis(i, p, u, U, d)
% INPUTS:
% - i : knot span index (Matlab index)
% - p : degree of the  B-spline
% - u : evaluation point
% - U : (u0, ..., um) clamped knot vector
% - d : maximum order of differentiation 
% ---------------------------------------------
% OUTPUT:
% - DN : Matrix (d+1) x (p+1).
%        The (k+1)-th row contains the derivatives of order k.
%        DN(1, :) -> Basis functions (order 0)
%        DN(2, :) -> First derivatives 
%        ...
% ---------------------------------------------

U = U(:)';                 % Force U to be a row vector 
d = min(d,p);              % It does not make sense to compute a derivative order higher than p 

% Pre-allocate
M  = zeros(p+1, p+1);       % M will contains all the basis fcts from degree 0 to degree p (on the upper triangular part)
                            % while in the lower triang. part it will contains the distances between two knots
DN = zeros(d+1, p+1);
left  = zeros(1,p);
right = zeros(1,p);

% M support matrix to store basis functions and knot
M(1, 1) = 1; % M(0,0) = 1

for j = 1:p
    % Assuming that i is the index of Matlab (given by the findspan
    % function) => u in [U(i), U(i+1))
    % U(i+1-j) corresponds to u_{i+1-j} and U(i+j) to u_{i+j} 
    left(j) = u - U(i+1-j);         
    right(j) = U(i+j) - u;         
    
    saved = 0;
    
    
    for r = 0:(j-1)
        % Lower triangle
        % Indeces of M start from 1
        M(j+1, r+1) = right(r+1) + left(j-r);   % Compute the denominator
        
        % Compute the term: N_i,j-1/distance
        if M(j+1, r+1)==0
            temp = 0;
        else
            temp = M(r+1, j) / M(j+1, r+1);    
        end
        
        % Upper triangle: compute the element N_i,j as the sum of the
        % previous 2 element N_i,j-1 and N_i-1,j-1
        M(r+1, j+1) = saved + right(r+1) * temp;    
        saved = left(j-r) * temp;                
    end
    
    M(j+1, j+1) = saved; 
end

% Load the basis functions (derivatives of order 0) in the first row of DN
for j = 0:p
    DN(1, j+1) = M(j+1, p+1); 
end

% This section computes the derivatives
for r = 0:p % Loop over function index 
    % s1 and s2 contains the 2 deriv. orders for computing the next one
    s1 = 1;  % Alternate rows in array a_k,j
    s2 = 2;
    
    a(1, 1) = 1; % a(0,0) = 1 
    
    % Loop to compute the k-th derivative
    for k = 1:d 
        dd = 0;      
        rk = r - k;
        pk = p - k;
        
        % The k-th derivative uses the basis elements from index r-k to r
        % => if r<k the index is 0 or negative => the elemnt basis =0
        if r >= k
            % Compute the coefficient a
            if M(pk+2, rk+1) == 0
                a(s2,1) = 0;
            else
                a(s2, 1) = a(s1, 1) / M(pk+2, rk+1); % M indices shifted +1
            end

            dd = a(s2, 1) * M(rk+1, pk+1);
        end

        if rk >= -1
            j1 = 1;
        else
            j1 = -rk;
        end
        
        if (r - 1) <= pk
            j2 = k - 1;
        else
            j2 = p - r;
        end
 
        for j = j1:j2 
            if M(pk+2, rk+j+1)==0
                a(s2,j+1) = 0;
            else
                % Note: j è indice logico (0..k), in Matlab usiamo j+1 per array 'a'
                a(s2, j+1) = (a(s1, j+1) - a(s1, j)) / M(pk+2, rk+j+1);
            end
            dd = dd + a(s2, j+1) * M(rk+j+1, pk+1);          % Update the sum
        end
        
        
        % The index of the basis fct must not exeed the number of available
        % basis fcts
        if r <= pk
            if M(pk+2, r+1) == 0
                a(s2,k+1)=0;
            else
                a(s2, k+1) = -a(s1, k) / M(pk+2, r+1);
            end           
            dd = dd + a(s2, k+1) * M(r+1, pk+1);            % Update the sum
        end
        
        DN(k+1, r+1) = dd;
        
        temp = s1;
        s1 = s2;
        s2 = temp;
    end
end

% Multiply through by the correct factors
r = p;
for k = 1:d
    for j = 0:p
        DN(k+1, j+1) = DN(k+1, j+1) * r;
    end
    r = r * (p - k);
end

end
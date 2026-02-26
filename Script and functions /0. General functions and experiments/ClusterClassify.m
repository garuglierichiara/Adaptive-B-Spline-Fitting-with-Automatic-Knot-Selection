%% THE FUNCTION ClusterClassify GROUPS THE ACTIVE KNOTS INTO CLUSTERS
% 
% It follows the algorithm 5: Cluster Classify from Kang et al. (2015)

function G = ClusterClassify(delta, U_tilde)
% INPUTS:
% - delta: threshold for constucting the clusters
% - U_tilde: knot vector (active knots)
% ---------------------------------------------
% OUTPUT:
% - G: vector containing group start indexes and end indexes
% ---------------------------------------------

k = 1;    % Initialize the cluster counter
G(1,1) = 1; % First index

% Loop over the element of U_tilde
for i = 2:numel(U_tilde)-1
    if U_tilde(i+1)-U_tilde(i)> delta   % if U(i+1) belongs to the next cluster
        G(k,2) = i;         % Save the end index of cluster k
        G(k+1,1)=i+1;       % Save the start index of the next cluster k+1
        k = k+1;            % Update the number of clusters
    end
end

G(k,2) = numel(U_tilde);    % Save the end index of the last cluster

G = reshape(G',2*k,1)';     % Reshape G from a matrix (kx2) to a vector (2kx1)
end

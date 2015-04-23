%% CrossQuery-Basic
function TopKResults = CQ_Basic(W, q, s, d, k, tilde_c, A_ID)

%%% Input parameters
%
% W: the transition matrix
% q: the ID of the query node of interest
% s: the ID of the source domain-specific network
% d: the ID of the target domain-specific network
% k: the number of retrieved nodes
% tilde_c: the normalized parameter for query preference
% A_ID: the IDs of domain nodes in each domain-specific network

%% Initialization

g = length(A_ID); % The number of domain-specific networks
DomainSizes = cellfun(@length,A_ID);

e = []; % Set initial scores

for i = 1:g
    
    tmp_e = zeros(DomainSizes(i),1);
    
    if i == s
        tmp_e(A_ID{i} == q) = 1; % Use g to replace 1 if NoN size is large
    end
    
    e = [e; tmp_e];
    
end

Wmax = max(W, [], 2);
Wmax = full(Wmax);

S = sum(DomainSizes(1:d-1))+1:sum(DomainSizes(1:d)); % The selected node set
S = S';

Iter = 1; % Iteration numer
p = e; % The random walk vector
Lower = (1-tilde_c)*p(S); % The lower bound vector

%% Score upper and lower bounds update loop
while length(S) > k
    
    % Update random walk scores
    %
    % Since MATLAB sparse matrix multiplication time complexity is
    % proportional to the number of non-zero entries, it is more efficient
    % in practice to directly multiply W and p than applying BFS search
    % before this multiplication in each iteration
    %
    % One can apply BFS one layer search in each iteration by using the
    % function BFS_layer
    p = W*p;
    
    % Update upper and lower bounds
    Lower = Lower + (1-tilde_c)*(tilde_c^Iter)*p(S);
    Upper = Lower + (tilde_c^(Iter+1))*Wmax(S);
    
    % Update threshold by the kth lower bound
    Theta = -kthvalue(-Lower, k);
    
    S = S(Upper >= Theta);
    Lower = Lower(Upper >= Theta);
    Upper = Upper(Upper >= Theta);
    
    % Avoid duplicates
    if max(Upper-Lower) < 1e-15
        SelectIdx = find(Lower > Theta);
        Duplicates = find(Lower == Theta);
        SelectIdx = [SelectIdx; Duplicates];
        SelectIdx = SelectIdx(1:k);
        S = S(SelectIdx);
    end
    
    Iter = Iter+1;
    
end

TopKResults = S;
TopKResults = TopKResults - sum(DomainSizes(1:d-1));
TopKResults = A_ID{d}(TopKResults);

end
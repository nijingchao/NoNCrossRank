%% CrossQuery-Fast
function [TopKResults, SubG_Idx] = CQ_Fast(Anorm, Y, G, q, s, d, k, alpha, c, epsilon, A_ID)

%%% Input parameters
%
% Anorm: the aggregated normalized adjacency matrix of domain-specific netowrks
% Y: the matrix encoding the cross-domain mapping information
% G: the adjacency matrix of the main network
% q: the ID of the query node of interest
% s: the ID of the source domain-specific network
% d: the ID of the target domain-specific network
% k: the number of retrieved nodes
% alpha: a regularization parameter for cross-network consistency
% c: a regularization parameter for query preference
% epsilon: an error factor to control the accuracy of results
% A_ID: the IDs of domain nodes in each domain-specific network

%% Initialization
DomainSizes = cellfun(@length,A_ID);
CumDomainSizes = cumsum(DomainSizes);

%% Extract relevant subnetwork from the main network
SubG_Idx = ExtractSubNet(G, s, d, epsilon);

%% Calculate matrices
g = length(SubG_Idx);
SubA_Idx = [];

for i = 1:g
    
    Hd = CumDomainSizes(SubG_Idx(i)) - DomainSizes(SubG_Idx(i)) + 1;
    Tl = CumDomainSizes(SubG_Idx(i));
    SubA_Idx = [SubA_Idx, Hd:Tl];
    
end

DomainSizes = DomainSizes(SubG_Idx);
A_ID = A_ID(SubG_Idx);
Anorm = Anorm(SubA_Idx',SubA_Idx);
s = find(SubG_Idx == s);
d = find(SubG_Idx == d);

SubG = G(SubG_Idx,SubG_Idx');
dSubG = sum(SubG,2);
Dy = cell(g,1);
for i = 1:g
    Dy{i} = dSubG(i)*speye(DomainSizes(i));
end
Dy = blkdiag(Dy{:});
Y = Y(SubA_Idx',SubA_Idx);
Dt = Dy - diag(sum(Y,2));
Y = Y + Dt;
Dyn = diag(diag(Dy).^(-0.5));
Ynorm = Dyn*Y*Dyn;

tilde_c = (c+2*alpha)/(1+2*alpha);
W = (c/(c+2*alpha))*Anorm + ((2*alpha)/(c+2*alpha))*Ynorm;

%% Apply CQ_Basic on the extracted NoN
TopKResults = CQ_Basic(W, q, s, d, k, tilde_c, A_ID);

end
%% Extract a relevant subnetwork from the main network w.r.t. source and target domains
function SubG_Idx = ExtractSubNet(G, s, d, epsilon)

%%% Input parameters
%
% G: the adjacency matrix of the main network
% s: the index of the source domain-specific network
% d: the index of the target domain-specific network
% epsilon: an error factor to control the accuracy of results

%% Intialiation for expansion
L_sd = Inf; % The shortest distance between s and d in G
MaxRadius = (L_sd-log10(epsilon))/2; % The maximal radius to search
Ns = []; % The neighbourhoods of s
Nd = []; % The neighbourhoods of d
rs = 0; % The radius of s
rd = 0; % The radius of d
SubG_Idx = []; % The indices of nodes in the extracted main network
Flag = 1; % The flag for the first overlap between neighbourhoods of s and d

%% Transform similarities in the main network to distances
dG = sum(G,2);
D_Gn = diag(dG.^(-0.5));
Gnorm = D_Gn*G*D_Gn;
Gnorm = max(Gnorm,eps);
DisG = -log10(Gnorm);
DisG = DisG - diag(diag(DisG));

%% Heap initialization
[rp, ci, vi] = sparse_to_csr(DisG);
n = length(rp)-1; % The number of nodes

Dis_s = Inf*ones(n,1); % The distances from s to other nodes
Hs = zeros(n,1); % The heap of node idices for s
Ps = zeros(n,1); % The heap positions of nodes for s
Len_s = 1; % The heap length of s
Hs(Len_s) = s;
Ps(s) = Len_s;
Dis_s(s) = 0;

Dis_d = Inf*ones(n,1); % The distances from d to other nodes
Hd = zeros(n,1); % The heap of node indices for d
Pd = zeros(n,1); % The heap positions of nodes for d
Len_d = 1; % The heap length of d
Hd(Len_d) = d;
Pd(d) = Len_d;
Dis_d(d) = 0;

%% Neighbourhood expansion loop
while (rs <= MaxRadius && Len_s > 0) || (rd <= MaxRadius && Len_d > 0)
    
    if length(intersect(Ns,Nd)) == 1 && Flag == 1 % The first overlap of two neighbourhoods
        
        Mid = intersect(Ns,Nd);
        L_sd = Dis_s(Mid)+Dis_d(Mid);
        MaxRadius = (L_sd-log10(epsilon))/2;
        Flag = 0;
        
    end
    
    if rs <= MaxRadius && Len_s > 0
        
        [NextNgbr_s, Len_s, Hs, Ps, Dis_s] = DijkstraExpansion(rp, ci, vi, Hs, Ps, Dis_s, Len_s);
        Ns = [Ns; NextNgbr_s];
        rs = Dis_s(NextNgbr_s);
        
    end
    
    if rd <= MaxRadius && Len_d > 0
        
        [NextNgbr_d, Len_d, Hd, Pd, Dis_d] = DijkstraExpansion(rp, ci, vi, Hd, Pd, Dis_d, Len_d);
        Nd = [Nd; NextNgbr_d];
        rd = Dis_d(NextNgbr_d);
        
    end
    
end

%% Full relax
% When the above loop stops, the shortest distances between s and the nodes
% in Nd but not in Ns are not fully computed. This also applies to the shortest
% distances between d and the nodes in Ns but not in Nd. Thus the following
% codes further relax these nodes.

Ns_extra = setdiff(Nd,Ns); % The nodes in Nd but not in Ns
Nd_extra = setdiff(Ns,Nd); % The nodes in Ns but not in Nd

while ~isempty(Ns_extra)
    
    [NextNgbr_s, Len_s, Hs, Ps, Dis_s] = DijkstraExpansion(rp, ci, vi, Hs, Ps, Dis_s, Len_s);
    if ismember(NextNgbr_s,Ns_extra)
        Ns = [Ns; NextNgbr_s];
        Ns_extra(Ns_extra==NextNgbr_s) = [];
    end
    
end

while ~isempty(Nd_extra)
    
    [NextNgbr_d, Len_d, Hd, Pd, Dis_d] = DijkstraExpansion(rp, ci, vi, Hd, Pd, Dis_d, Len_d);
    if ismember(NextNgbr_d,Nd_extra)
        Nd = [Nd; NextNgbr_d];
        Nd_extra(Nd_extra==NextNgbr_d) = [];
    end
    
end

%% Further prunning Ns and Nd
Nsd = [Ns; Nd]; % The combined neighbourhood of s and d
Nsd = unique(Nsd);

for i = 1:length(Nsd)
    
    u = Nsd(i);
    if Dis_s(u) + Dis_d(u) <= L_sd - log10(epsilon)
        SubG_Idx = [SubG_Idx; u];
    end
    
end

end
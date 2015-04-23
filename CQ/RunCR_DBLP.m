%% CrossRank evaluation on DBLP dataset
function TopKAuthorNames = RunCR_DBLP(alpha, c, MaxIter, epsilon, q, s, d, k)

%%% Input parameters
%
% If no input parameters are provided, the default values will be used.
%
% alpha: a regularization parameter for cross-network consistency
% c: a regularization parameter for query preference
% MaxIter: the maximal number of iteration for updating ranking vector
% epsilon: a convergence parameter
% q: the ID of the query node of interest
% s: the ID of the source domain-specific network
% d: the ID of the target domain-specific network
% k: the number of retrieved nodes
%
% Looking at ConfDict in ../ExampleDatasets/DBLP_NoN.mat to determine
% source and target domain IDs, s and d.
%
% Looking at AuthorDict in ../ExampleDatasets/DBLP_NoN.mat to determine the
% ID of the query node q.
%
% RunCR_DBLP returns the names of top k relevant authors from the target
% domain-specific network

%% Parameter initialization
if nargin < 8
    k = 10;
end
if nargin < 7
    d = 20; % The ID of SIGMOD Conference
end
if nargin < 6
    s = 1; % The ID of KDD
end
if nargin < 5
    q = 121; % The ID of Jiawei Han
end
if nargin < 4
    epsilon = 1e-15;
end
if nargin < 3
    MaxIter = 1000;
end
if nargin < 2
    c = 0.85;
end
if nargin < 1
    alpha = 0.2;
end

%% Load NoN data
load('../ExampleDatasets/DBLP_NoN.mat');

%% Rename networks
g = length(CoAuthorNets); % The number of domain-specific networks
G = ConfNet; % The main network
A = CoAuthorNets; % The domain-specific networks
A_ID = CoAuthorNetsID; % The IDs of nodes in domain-specific networks

%% CR Precomputation, this step only needs to be done once for a dataset
PrecompFileName = 'Precomp_Values_DBLP.mat';

if exist(PrecompFileName, 'file') == 2
    disp('A precomputation file has been detected ...');
else
    disp('Precomputation starts ...');
    Precomputation(A, A_ID, G, PrecompFileName);
end

disp('Load the precomputation file ...');
load(PrecompFileName);

%% Run CR
% Set initial scores
tic;
DomainSizes = cellfun(@length,A_ID);
e = [];

for i = 1:g
    
    tmp_e = zeros(DomainSizes(i),1);
    
    if i == s
        tmp_e(A_ID{i} == q) = 1; % Use g to replace 1 if NoN size is large
    end
    
    e = [e; tmp_e];
    
end

% CR
[r, Objs, Deltas] = CR(Anorm, Ynorm, I_n, e, alpha, c, MaxIter, epsilon);

% Sort ranking scores in the target domain-specific network
rd_Idx = sum(DomainSizes(1:d-1))+1:sum(DomainSizes(1:d));
rd = r(rd_Idx);
[Sort_rd, Sort_Idx] = sort(rd,'descend');
TopKResults = Sort_Idx(1:k);
TopKResults = A_ID{d}(TopKResults);
RunTime = toc;

disp(['The running time of CR is ' num2str(RunTime) ' seconds.']);

TopKAuthorNames = AuthorDict(TopKResults);

end
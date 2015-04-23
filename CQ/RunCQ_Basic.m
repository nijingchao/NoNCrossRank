%% CrossQuery-Basic evaluation on DBLP dataset
function TopKAuthorNames = RunCQ_Basic(alpha, c, q, s, d, k)

%%% Input parameters
%
% If no input parameters are provided, the default values will be used.
%
% alpha: a regularization parameter for cross-network consistency
% c: a regularization parameter for query preference
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
% RunCQ_Basic returns the names of top k relevant authors from the target
% domain-specific network

%% Parameter initialization
if nargin < 6
    k = 10;
end
if nargin < 5
    d = 20; % The ID of SIGMOD Conference
end
if nargin < 4
    s = 1; % The ID of KDD
end
if nargin < 3
    q = 121; % The ID of Jiawei Han
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
G = ConfNet; % The main network
A = CoAuthorNets; % The domain-specific networks
A_ID = CoAuthorNetsID; % The IDs of nodes in domain-specific networks

%% CQ Precomputation, this step only needs to be done once for a dataset
PrecompFileName = 'Precomp_Values_DBLP.mat';

if exist(PrecompFileName, 'file') == 2
    disp('A precomputation file has been detected ...');
else
    disp('Precomputation starts ...');
    Precomputation(A, A_ID, G, PrecompFileName);
end

disp('Load the precomputation file ...');
load(PrecompFileName);

%% Run CQ_Basic
% Initialization
tic;
tilde_c = (c+2*alpha)/(1+2*alpha);
W = (c/(c+2*alpha))*Anorm + ((2*alpha)/(c+2*alpha))*Ynorm;

% CQ_Basic
TopKResults = CQ_Basic(W, q, s, d, k, tilde_c, A_ID);
RunTime = toc;

disp(['The running time of CQ_Basic is ' num2str(RunTime) ' seconds.']);

TopKAuthorNames = AuthorDict(TopKResults);

end
%% CrossQuery-Fast evaluation on DBLP dataset
function [TopKAuthorNames, RelevantDomains] = RunCQ_Fast(alpha, c, epsilon, q, s, d, k)

%%% Input parameters
%
% If no input parameters are provided, the default values will be used.
%
% alpha: a regularization parameter for cross-network consistency
% c: a regularization parameter for query preference
% epsilon: an error factor to control the accuracy of results
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
% RunCQ_Fast returns the names of top k relevant authors from the target
% domain-specific network and the relevant domains of the source and target
% domains

%% Parameter initialization
if nargin < 7
    k = 10;
end
if nargin < 6
    d = 20; % The ID of SIGMOD Conference
end
if nargin < 5
    s = 1; % The ID of KDD
end
if nargin < 4
    q = 121; % The ID of Jiawei Han
end
if nargin < 3
    epsilon = 0.003;
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

%% Run CQ_Fast
% Initialization
tic;
% CQ_Fast
[TopKResults, SubG_Idx] = CQ_Fast(Anorm, Y, G, q, s, d, k, alpha, c, epsilon, A_ID);
RunTime = toc;

disp(['The running time of CQ_Fast is ' num2str(RunTime) ' seconds.']);

TopKAuthorNames = AuthorDict(TopKResults);
RelevantDomains = ConfDict(SubG_Idx);

end
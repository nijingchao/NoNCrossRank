%% CR leave-one-out cross validation on tissue-specific PPI networks
function CR_CrossValidation(alpha, c, MaxIter, epsilon)

%%% Input parameters
%
% If no input parameters are provided, the default values will be used.
%
% alpha: a regularization parameter for cross-network consistency
% c: a regularization parameter for query preference
% MaxIter: the maximal number of iteration for updating ranking vector
% epsilon: a convergence parameter

%% Parameter initialization
if nargin < 4
    epsilon = 1e-6;
end
if nargin < 3
    MaxIter = 1000;
end
if nargin < 2
    c = 0.85;
end
if nargin < 1
    alpha = 0.5;
end

%% Load NoN data
load('../ExampleDatasets/P_G_NoN.mat');

%% Rename networks
g = length(TSGeneNets); % The number of domain-specific networks
G = PhenotypeSimNet; % The main network
A = TSGeneNets; % The domain-specific networks
A_ID = TSGeneNetsID; % The IDs of nodes in domain-specific networks
A_Seeds = Seeds; % The seed/query nodes in domain-specific networks

%% CR precomputation, this step only needs to be done once for a dataset
PrecompFileName = 'CR_Precomp_Values_TPPI.mat';

if exist(PrecompFileName, 'file') == 2
    disp('A precomputation file has been detected ...');
else
    disp('CR precomputation starts ...');
    CR_Precomputation(A, A_ID, G, PrecompFileName);
end

disp('Load the precomputation file ...');
load(PrecompFileName);

%% Leave-one-out cross validation
% Expand test genes s.t. test gene one by one
ExpandSeeds = vertcat(A_Seeds{:});

% Leave-one-out cross validation loop
RankScoreRecord = cell(1,length(ExpandSeeds));
RankRecord = cell(1,length(ExpandSeeds));
TotalCounter = 0;

disp('Leave-one-out cross validation starts ...');

for j = 1:g
    
    for t = 1:length(A_Seeds{j})
        
        TotalCounter = TotalCounter + 1;
        
        % Initialize query vector
        e = [];
        
        for i = 1:g
            
            if i == j
                head = length(e) + 1; % Record head position in e/r for final evaluation
            end
            
            tmp_e = zeros(length(A_ID{i}),1);

            if ismember(A_Seeds{j}(t),A_Seeds{i})
                subn = length(A_Seeds{i}) - 1;
                if subn ~= 0
                    [Fia, seedidx] = ismember(A_Seeds{i},A_ID{i});
                    tmp_e(seedidx) = 1/subn;
                    [Fia1,seedidx1] = ismember(A_Seeds{j}(t),A_ID{i});
                    tmp_e(seedidx1) = 0;
                end
            else
                subn = length(A_Seeds{i});
                [Fia, seedidx] = ismember(A_Seeds{i},A_ID{i});
                tmp_e(seedidx) = 1/subn;
            end
            
            e = [e; tmp_e];
            
            if i == j
                tail = length(e); % Record tail position in e/r for final evaluation
            end
            
        end
        
        % CR
        e = sparse(e);
        [r, Objs, Deltas] = CR(Anorm, Ynorm, I_n, e, alpha, c, MaxIter, epsilon);
        
        % Record results
        RankScore = r(head:tail);
        RankScore = (round(RankScore*1e16))/1e16;
        FullRankScore = zeros(length(AllGeneID),1);
        [Fia, idx] = ismember(A_ID{j}, AllGeneID);
        FullRankScore(idx) = RankScore;
        
        SeedGeneList = setdiff(A_Seeds{j}, A_Seeds{j}(t));
        [proj, ia, idx] = intersect(SeedGeneList, AllGeneID);
        FullRankScore(idx) = 0;
        
        RankScoreRecord{TotalCounter} = FullRankScore;
        
        % Sort results
        [SortFullRankScore, IX] = sort(FullRankScore, 'descend');
        RankRecord{TotalCounter} = IX;
        
        disp(['Finished Number of Folds/Total Number of Folds: ' num2str(TotalCounter) '/' num2str(length(ExpandSeeds))]);
        
    end
    
end

%% Save results
disp('Leave-one-out cross validation finishes, save results ...');

save('CRResults.mat','RankScoreRecord','RankRecord','ExpandSeeds','AllGeneID');

%% Evaluation
disp('AUC value evaluation starts ...');

AUCEvaluation(RankRecord, ExpandSeeds, AllGeneID);

end
%% CR precomputation
function CR_Precomputation(A, A_ID, G, PrecompFileName)

%%% Input parameters
%
% A: the domain-specific networks
% A_ID: the corresponding IDs of domain-specific networks in A
% G: the adjacency matrix of the main network
% PrecompFileName: the file name to store precomputation results

%% Initialization
g = length(A);
ns = cellfun(@length, A_ID);
n = sum(ns);
I_n = speye(n);

%% Normalize A
Anorms = cell(size(A));

for i = 1:g
    
    D = sum(A{i},2);
    D = D.^(-0.5);
    D = diag(D);
    Anorms{i} = D*A{i}*D;
    
end

Anorm = blkdiag(Anorms{:});

%% Common node mapping
Y = sparse(n,n);
Dy = cell(g,1);
dG = sum(G,2);

for i = 1:g
    
    Dy{i} = dG(i)*speye(ns(i));
    Y_i = sparse(ns(i),n);
    
    for j = i:g
        
        [proj, I1, I2] = intersect(A_ID{i}, A_ID{j});
        Oij = sparse(I1, I2, ones(length(proj),1), ns(i), ns(j));
        Y_i(:, sum(ns(1:j-1))+1:sum(ns(1:j))) = G(i,j)*Oij;
        
    end
    
    Y(sum(ns(1:i-1))+1:sum(ns(1:i)), :) = Y_i;
    
end

Y = triu(Y) + triu(Y)' - diag(diag(Y));

%% Construct normalized Y
Dy = blkdiag(Dy{:});
Dt = Dy - diag(sum(Y,2));
Y = Y + Dt;
Dyn = diag(diag(Dy).^(-0.5));
Ynorm = Dyn*Y*Dyn;

%% Save precomputed matrices
save(PrecompFileName, 'Anorm', 'Ynorm', 'I_n');

end
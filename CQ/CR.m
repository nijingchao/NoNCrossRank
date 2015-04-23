%% CR
function [r, Objs, Deltas] = CR(Anorm, Ynorm, I_n, e, alpha, c, MaxIter, epsilon)

%%% Input parameters
%
% Anorm: the aggregated normalized adjacency matrix of domain-specific netowrks
% Ynorm: the normalized matrix encoding the cross-domain mapping information
% I_n: an identity matrix of size n
% e: the query vector
% c: a regularization parameter for query preference.
% MaxIter: the maximal number of iteration for updating ranking vector.
% epsilon: a convergence parameter.

%% Initialization
% Initialize ranking vector
r = e;

% Initialize parameters
gamma = c/(1 + 2*alpha);
kappa = 2*alpha/(1 + 2*alpha);
eta = (1 - c)/(1 + 2*alpha);

% Convergence analysis parameters
% Either the difference between objective values or ranking vector norms
% can be used as a measure of convergence.

% J1 = J_CR(Anorm, Ynorm, I_n, r, e, alpha, c); % Objective value measure
J1 = r; % ranking vector norm measure
J1 = (round(J1*1e16))/1e16;
delta = 99999;
Objs = [];
Deltas = [];
Iter = 1;

%% Power method update loop
while delta > epsilon && Iter <= MaxIter
    
    M = gamma*Anorm+kappa*Ynorm;
    r = M*r+eta*e;
    
    % Convergence analysis
    J2 = J1;
%     J1 = J_CR(Anorm, Ynorm, I_n, r, e, alpha, c); % Objective value measure
    J1 = r; % ranking vector norm measure
    J1 = (round(J1*1e16))/1e16;
%     delta = J2 - J1; % Objective value measure
    delta = norm(J2 - J1, 1); % ranking vector norm measure
    
%     Objs = [Objs, J1]; % Objective value measure
    Objs = [Objs, norm(J1, 1)]; % ranking vector norm measure
    Deltas = [Deltas, delta];
    
    Iter = Iter + 1;
    
end

end
%% CrossRank objective function value
function Obj = J_CR(Anorm, Ynorm, I_n, r, e, alpha, c)

X = I_n - Ynorm;
Obj = r'*c*(I_n-Anorm)*r + (1-c)*norm(r-e, 'fro')^2 + 2*alpha*(r'*X*r);

end
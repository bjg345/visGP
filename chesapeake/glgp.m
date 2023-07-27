
grid_train = readmatrix('./train_loc.csv');
grid_test = readmatrix('./test_loc.csv');
y = readmatrix('./train_y.csv');
vals_test = readmatrix('./test_y.csv');

X1 = grid_train;
X = grid_test;

L_cur = inf
param_cur = NaN;
for k = 50:50:150
    for eps = exp(-2:.2:2)
        fun = @(par) min(-1*CovMatrix_likelihood(X1, X, y, k, eps, par(1), par(2)), 9e99);
        [param, L]= fmincon(fun, [1,1], [], [], [], [], [0,0], []);
        if L < L_cur
            k_cur = k;
            eps_cur = eps;
            L_cur = L;
            param_cur = param;
        end
    end
end

param = [k_cur, eps_cur, param_cur];

vals_pred = CovMatrix_predict(X1, X, y, param(1), param(2), param(3), param(4));

save('glgp_pred.mat', 'vals_pred')


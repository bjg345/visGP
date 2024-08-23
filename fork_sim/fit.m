function [fit] = fit(data1,data2,y)


%   Inputs (required)
%       data1   : data matrix for the labeled points 
%       data2   : data matrix for the unlabeled points
%       y  :  observed values of the regression function f over data1

%   Outputs
%     optimal parameter values

data=[data1;data2];

[n, ~] = size(data) ;
[n1, ~] = size(data1) ;
[n2, ~] = size(data2) ;

eps_space = [.001, .0025, .005, .0075, .01, .025, .05, .075, .1, .25, .5, .75, 1, 2.5, 5]  ;
sig_space = eps_space
k_space = [50];
t_space = [5:5:35];


cur_val = -inf;
cur_k = nan;
cur_eps = nan;
cur_t = nan;
cur_sig = nan;
for eps = eps_space
    for k = k_space
        for t = t_space
            for sig = sig_space
		  [k, eps, t, sig]
                  val = CovMatrix_likelihood(data1, data2, y, k, eps, ...
                      t, sig)
                 if val>cur_val
                     cur_val = val;
                     cur_k = k;
                     cur_eps = eps;
                     cur_t = t;
                     cur_sig = sig;
                 end
            end
         end
    end
end

fit = [cur_k, cur_eps, cur_t, cur_sig, cur_val]

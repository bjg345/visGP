function [L] = CovMatrix_likelihood(data1,data2,y,k,eps,t,sig)


%   Inputs (required)
%       data1   : data matrix for the labeled points 
%       data2   : data matrix for the unlabeled points
%       y  :  observed values of the regression function f over data1
%       k    : first k eigenpairs of the Graph Laplacian over data1 and
%       data2
%       eps  : bandwidth in the Graph Laplacian
%       t  :  diffusion time
%       sig : covariance of the the error
%   Outputs
%      L : the loglikelihood of fit

data=[data1;data2];

[n, ~] = size(data) ;
[n1, ~] = size(data1) ;
[n2, ~] = size(data2) ;

[index,distance]= knnsearch(data, data,'k', n);
distance(:,1)=0;
ker = exp(-distance.^2/(4*eps));
ii = (1:n)'*ones(1,n);
W = sparse(ii, index, ker, n, n);
D = sum(W, 2);
W = bsxfun(@rdivide, bsxfun(@rdivide, W, D), transpose(D));
D = sqrt(sum(W, 2)); 
W = bsxfun(@rdivide, bsxfun(@rdivide, W, D), transpose(D));
[U,S] = eigs(W, k);
if(isnan(S(1,1)) | ~isreal(S(k,k)))
	L = -9e99;
	return
end
U = bsxfun(@rdivide, U(:, 1:end), D);
% Construct heat kernel
C = zeros(n,n);
for i=1:k
C = C + exp(t*(S(i,i)-1)/eps)*U(1:n,i)*U(1:n,i)'/norm(U(1:n,i))^2*n;
end
% log-likelihood
try 
    L = log(mvnpdf(y-mean(y), [], C(1:n1, 1:n1)+sig*eye(n1)));
catch 
    L = -9e99;
end


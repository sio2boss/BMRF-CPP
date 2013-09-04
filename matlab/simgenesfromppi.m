% simulate gene expression data using GG model and MRF model
% reference: A Markov random field model for network-based analysis of
% genomic data by Zhi Wei and Hongze Li
% Input: ppi, X0, m, n
% output: gene expression data Y

function [X Y] = simgenesfromppi(ppi, X0, m, n, weight)

gid = unique(ppi(:));
% step 1, simulate X, perform sampling based on X0, according to MRF
% conditional likelihood
gamma0 = 10; gamma1 = 10; beta = 20;
delta = 0;
for i=1:length(gid)
    gconn = union(ppi(find(ppi(:,1)==gid(i)),2), ppi(find(ppi(:,2)==gid(i)),1));
    [a,b] = intersect(gid, gconn);
    u1 = (X0(i)*weight + sum(X0(b)==1))/(weight+length(gconn));    
    u0 = ((1-X0(i))*weight + sum(X0(b)==0))/(weight+length(gconn));
    UU1(i) = u1;
    UU0(i) = u0;
    mrfpdf(i) = exp((1-X0(i))*(gamma0-beta*u1)+X0(i)*(gamma1-beta*u0))/(exp(gamma0-beta*u1)+exp(gamma1-beta*u0));
end
% delta = 0.6;
for i=1:length(gid)
    if X0(i) == 1
        p1 = max(1-mrfpdf(i)-delta, 0);
        p2 = min(mrfpdf(i) + delta, 1);
        X(i) = randsample([0 1],1,true,[p1 p2]);
    else
        p1 = min(mrfpdf(i) + delta, 1);
        p2 = max(1-mrfpdf(i)-delta, 0);
        X(i) = randsample([0 1], 1, true, [p1 p2]);
    end
end

% step 2, simulate gene expression level Y according to GG model
alpha = 10; alpha0 = 0.9; v = 0.5;
% K1 = v^alpha0*gamma(m*alpha+alpha0)/(gamma(alpha)^m*gamma(alph0));
% K2 = v^alpha0*gamma(n*alpha+alpha0)/(gamma(alpha)^n*gamma(alph0));
% K = v^alpha0*gamma((m+n)*alpha+alpha0)/(gamma(alpha)^(m+n)*gamma(alph0));
index0 = find(X == 0);
lambda = gamrnd(alpha0, 1./v, length(index0), 1);
for i=1:(m+n)
    Y(index0,i) = gamrnd(alpha, 1./lambda);
end
index1 = find(X == 1);
lambda1 = gamrnd(alpha0, 1./v, length(index1),1);
lambda2 = gamrnd(alpha0, 1./v, length(index1), 1);
for i = 1:m
    Y(index1,i) = gamrnd(alpha, 1./lambda1);
end
for i = (m+1):(m+n)
    Y(index1, i) = gamrnd(alpha, 1./lambda2);
end

    
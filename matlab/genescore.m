function [ Zscore ] = genescore( x,y )
%TTEST2 Two-sample t-test with pooled or unpooled variance estimate.
%   H = TTEST2(X,Y) performs a t-test of the hypothesis that two
%   independent samples, in the vectors X and Y, come from distributions
%   with equal means, and returns the result of the test in H.  H=0
%   indicates that the null hypothesis ("means are equal") cannot be
%   rejected at the 5% significance level.  H=1 indicates that the null
%   hypothesis can be rejected at the 5% level.  The data are assumed to
%   come from normal distributions with unknown, but equal, variances.  X
%   and Y can have different lengths.
%
%   This function performs an unpaired two-sample t-test. For a paired
%   test, use the TTEST function.
%
%   X and Y can also be matrices or N-D arrays.  For matrices, TTEST2
%   performs separate t-tests along each column, and returns a vector of
%   results.  X and Y must have the same number of columns.  For N-D
%   arrays, TTEST2 works along the first non-singleton dimension.  X and Y
%   must have the same size along all the remaining dimensions.
%
%   TTEST2 treats NaNs as missing values, and ignores them.
%
%   H = TTEST2(X,Y,ALPHA) performs the test at the significance level
%   (100*ALPHA)%.  ALPHA must be a scalar.
%
%   H = TTEST2(X,Y,ALPHA,TAIL) performs the test against the alternative
%   hypothesis specified by TAIL:
%       'both'  -- "means are not equal" (two-tailed test)
%       'right' -- "mean of X is greater than mean of Y" (right-tailed test)
%       'left'  -- "mean of X is less than mean of Y" (left-tailed test)
%   TAIL must be a single string.
%
%   H = TTEST2(X,Y,ALPHA,TAIL,VARTYPE) allows you to specify the type of
%   test.  When VARTYPE is 'equal', TTEST2 performs the default test
%   assuming equal variances.  When VARTYPE is 'unequal', TTEST2 performs
%   the test assuming that the two samples come from normal distributions
%   with unknown and unequal variances.  This is known as the Behrens-Fisher
%   problem. TTEST2 uses Satterthwaite's approximation for the effective
%   degrees of freedom.  VARTYPE must be a single string.
%
%   [H,P] = TTEST2(...) returns the p-value, i.e., the probability of
%   observing the given result, or one more extreme, by chance if the null
%   hypothesis is true.  Small values of P cast doubt on the validity of
%   the null hypothesis.
%
%   [H,P,CI] = TTEST2(...) returns a 100*(1-ALPHA)% confidence interval for
%   the true difference of population means.
%
%   [H,P,CI,STATS] = TTEST2(...) returns a structure with the following fields:
%      'tstat' -- the value of the test statistic
%      'df'    -- the degrees of freedom of the test
%      'sd'    -- the pooled estimate of the population standard deviation
%                 (for the equal variance case) or a vector containing the
%                 unpooled estimates of the population standard deviations
%                 (for the unequal variance case)
%
%   [...] = TTEST2(X,Y,ALPHA,TAIL,VARTYPE,DIM) works along dimension DIM of
%   X and Y.  Pass in [] to use default values for ALPHA, TAIL, or VARTYPE.
%
%   See also TTEST, RANKSUM, VARTEST2, ANSARIBRADLEY.

%   References:
%      [1] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, section 13.4. (Table 13.4.1 on page 210)

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $  $Date: 2011/05/09 01:27:10 $

if nargin < 2
    error(message('stats:ttest2:TooFewInputs'));
end

if nargin < 6
    % Figure out which dimension mean will work along by looking at x.  y
    % will have be compatible. If x is a scalar, look at y.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = find(size(y) ~= 1, 1); end
    if isempty(dim), dim = 1; end
    
    % If we haven't been given an explicit dimension, and we have two
    % vectors, then make y the same orientation as x.
    if isvector(x) && isvector(y)
        if dim == 2
            y = y(:)';
        else % dim == 1
            y = y(:);
        end
    end
end

% Make sure all of x's and y's non-working dimensions are identical.
sizex = size(x); sizex(dim) = 1;
sizey = size(y); sizey(dim) = 1;
if ~isequal(sizex,sizey)
    error(message('stats:ttest2:InputSizeMismatch'));
end

xnans = isnan(x);
if any(xnans(:))
    nx = sum(~xnans,dim);
else
    nx = size(x,dim); % a scalar, => a scalar call to tinv
end
ynans = isnan(y);
if any(ynans(:))
    ny = sum(~ynans,dim);
else
    ny = size(y,dim); % a scalar, => a scalar call to tinv
end

s2x = nanvar(x,[],dim);
s2y = nanvar(y,[],dim);
difference = nanmean(x,dim) - nanmean(y,dim);
% equal variances
dfe = nx + ny - 2;
sPooled = sqrt(((nx-1) .* s2x + (ny-1) .* s2y) ./ dfe);
se = sPooled .* sqrt(1./nx + 1./ny);
ratio = difference ./ se;
test = -abs(ratio);

% Compute the correct p-value for the test, and confidence intervals
% if requested. two-tailed test
p = 2 * tcdf(test,dfe);
p = 1-p;

% Fixed version
machine_epsilon = 1e-16;
p(p==1.0) = 1.0 - machine_epsilon;
p(p<machine_epsilon) = machine_epsilon;

% Old version that induces boundary issues
% p(p==1.0) = 0.9999;
% p(p<1e-15) = 1e-15;

gene_Zscore = icdf('normal', p, 0,1);
Zscore = (gene_Zscore-mean(gene_Zscore))/std(gene_Zscore);

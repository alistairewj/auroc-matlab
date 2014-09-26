function [ Wvar, W ] = wilcoxonVariance(pred, target)
%WILCOXONVARIANCE	Calculates the variance of the Wilcoxon statistic given predictions and targets
%	[ Wvar, W ] = wilcoxonVariance(pred, target)
%
%
%	Inputs:
%       pred   - the prediction made (the probability that target=1).
%       target - the target variable, 0 or 1, predicted by pred
%
%
%
%	Outputs:
%       Wvar - Variance in the Wilcoxon statistic
%		W    - Probability( PRED|target=1 > PRED|target=0 )
%               Calculation: sum(sum(PRED|target=1 > PRED|target=0))
%               Equivalent to the AUROC.
%
%
%	Example
%		[ Wvar, W ] = wilcoxonVariance(pred, target)
%
%	See also FCN1

%	References:
%
%

%	Copyright 2011 Alistair Johnson

%   $LastChangedBy: alistair $
%   $LastChangedDate: 2014-07-30 10:31:51 +0100 (Wed, 30 Jul 2014) $
%   $Revision: 1182 $
%   Originally written on GLNXA64 by Alistair Johnson, 17-Oct-2011 15:37:51
%   Contact:alistairewj@gmail.com


% Split into positive and negative vectors
idx=target==1;
positive=sort(pred(idx),'ascend'); negative=sort(pred(~idx),'descend');
N1=length(positive); N0=length(negative);

% First calculate "placements", V0 and V1
V0=zeros(size(negative)); V1=zeros(size(positive));

[loopVar,loopVarIdx]=sort([N0,N1]); % efficient for loops-> no redundancy

for k=1:loopVar(1) % Initially loop over both positive/negative
    % Using the fact that positive/negative are sorted...
    
    % V1: Find proportion of negative each point in positive exceeds
    % V0: Find proportion of positive each point in negative is less than
end

if loopVarIdx(1)==1 % More negative outcomes than positives
    for k=loopVar(1)+1:loopVar(2) % Continue looping over negative
        
    end
elseif loopVarIdx(1)==2 % More positive outcomes than positives
    for k=loopVar(1)+1:loopVar(2) % Continue looping over positive
        
    end
end

Wvar=var(V0)/N0 + var(V1)/N1; % Approximation: since N0*N1 >> 10,000, W*(1-W)/(N0*N1) ~ 0


end
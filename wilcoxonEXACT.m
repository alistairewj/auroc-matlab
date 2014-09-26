function [ W ] = wilcoxonEXACT( pred, target )
%WILCOXON Calculates the Wilcoxon statistic, equivalent to the area under
%   the receiver operator characteristic (AUROC)
%
%   W = wilcoxon(pred, target) Calculates the Wilcoxon statistic W given a
%   vector of predictions ranging from 0-1 and their associates targets,
%   which are either 0 or 1
%
%	Inputs:
%       target - the target variable, 0 or 1, predicted by pred
%       pred   - the prediction made (the probability that target=1).
%   pred,target should be column vectors.
%
%	Outputs:
%		W - Probability( PRED|target=1 > PRED|target=0 )
%               Calculation: sum(sum(PRED|target=1 > PRED|target=0))
%               Equivalent to the AUROC.
%
%
%   Example 1: Classify ~90% correctly.
%       pred=rand(100,1);
%       target=round(pred);
%       target(1:10:end)=1-target(1:10:end);
%       W = wilcoxonEXACT(pred,target)

%	Copyright 2011 Alistair Johnson

%   $LastChangedBy: alistair $
%   $LastChangedDate: 2014-07-30 12:14:35 +0100 (Wed, 30 Jul 2014) $
%   $Revision: 1183 $
%   Originally written on GLNXA64 by Alistair Johnson, 17-Oct-2011 12:11:42
%   Contact:alistairewj@gmail.com

% do not use nan targets in calculations
idxNotNan = isfinite(target);
idx=target==0;

P = size(pred,2);

negative=pred(idx & idxNotNan,:);
N1=size(negative,1);
positive=pred(~idx & idxNotNan,:);
N2=size(positive,1);

W=zeros(1,P);
W2=zeros(1,P);
for p=1:P
    for n=1:N2
        % compare 0s to 1s
        idx=gt(positive(n,p),negative(:,p));
        W(p)=W(p)+sum(idx,1);
        
        % Requires double computation time, but more accurate.
        W2(p)=W2(p)+sum(eq(positive(n,p),negative(:,p))); 
    end
end
W=W+W2(p)*0.5;
W=W/(N1*N2);
end

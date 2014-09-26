function [ W ] = wilcoxonEXACTXY( X, Y )
%WILCOXON Calculates the Wilcoxon statistic, equivalent to the area under
%   the receiver operator characteristic (AUROC)
%
%   W = wilcoxon(X, Y) Calculates the Wilcoxon statistic W given a
%   vector of X predictions corresponding to a target value of 1 and a
%   vector of Y predictions corresponding to a target value of 0. If X and
%   Y are matrices, then the W statistic is calculated for each column of X
%   and Y and returned as a row vector.
%
%	Inputs:
%       X   - predictions for target value of 1
%       Y   - predictions for target value of 0
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
%       X=pred(target==1);
%       Y=pred(target==0);
%       W = wilcoxonEXACTXY(X,Y)

%	Copyright 2011 Alistair Johnson

%   $LastChangedBy: alistair $
%   $LastChangedDate: 2014-07-30 10:31:51 +0100 (Wed, 30 Jul 2014) $
%   $Revision: 1182 $
%   Originally written on GLNXA64 by Alistair Johnson, 17-Oct-2011 12:11:42
%   Contact:alistairewj@gmail.com

N1=size(Y);
N2=size(X);

if N1(2)~=N2(2) % Matrices do not have predictions for both target values
    error('X and Y should have the same number of columns.\n');
end

K=N1(2); % Number of X/Y pairs
W=zeros(1,K); W2=zeros(1,K);
for r=1:K % across columns, i.e. different predictors
    % compare 0s to 1s
    for n=1:N2(1)
        idx=gt(X(n,r),Y(:,r));
        W(r)=W(r)+sum(idx);
        % W2 requires double computation time, but more accurate.
        W2(r)=W2(r)+sum(eq(X(n,r),Y(:,r)));
    end
    W(r)=W(r)+W2(r)*0.5;
    W(r)=W(r)/(N1(1)*N2(1));
end
end
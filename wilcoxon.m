function [ W ] = wilcoxon( pred, target )
%WILCOXON Calculates the Wilcoxon statistic, equivalent to the area under
%   the receiver operator characteristic (AUROC)
%
%   W = wilcoxon(pred, target) Calculates the Wilcoxon statistic W given a
%   vector of predictions ranging from 0-1 and their associates targets,
%   which are either 0 or 1
%   
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
%       W = wilcoxon(pred,target)

%	$LastChangedBy: alistair $
%	$LastChangedDate: 2014-07-30 10:31:51 +0100 (Wed, 30 Jul 2014) $
%	$Revision: 1182 $
%	Originally written on GLNXA64 by Alistair Johnson, 17-Oct-2011 12:11:42
%	Contact: alistairewj@gmail.com

%=== Arrange predictions
[pred,idx] = sort(pred,1,'ascend');
target=target(idx);
[N,P] = size(pred);

%=== Find location of negative targets
W=zeros(1,P); % 1xP where P is # of AUROCs to calculate
negative=false(N,P);

for n=1:P
    negative(:,n) = target(:,n)==0;
    
    %=== Count the number of negative targets below each element
    negativeCS = cumsum(negative(:,n),1);
    
    %=== Only get positive targets
    pos = negativeCS(~negative(:,n));
    
    W(n) = sum(pos);
end
count=sum(negative,1); %=== count number who are negative
count = count .* (N-count); %=== multiply by positives
W=W./count;

end

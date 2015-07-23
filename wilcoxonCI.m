function [ auc, auc_ci, S, S10, S01, V10, V01 ] = wilcoxonCI(pred,target, alpha)
%WILCOXONCOVARIANCE	Calculates the covariance for R sets of predictions and outcomes
%	[ auc, auc_ci, S, S10, S01, V10, V01 ] = wilcoxonCI(pred,target, alpha)
%	
%
%	Inputs: Must come in pairs in the following order:
%       pred   - the prediction made (the probability that target=1).
%       target - the target variable, 0 or 1, predicted by pred
%       alpha  - Critical value for confidence intervals (default 0.05)
%
%	Outputs:
%		S       - Covariance matrix of the AUROCs
%		S10     - Covariance matrix for positive outcomes
%		S01     - Covariance matrix for negative outcomes
%		V10     - Pairwise comparisons of positive predictions
%		V01     - Pairwise comparisons of negative predictions
%       auc     - Row vector containing AUROC/Wilcoxon statistic
%		
%
%	Example
%		[ auc, auc_ci, S, S10, S01, V10, V01 ] = ...
%           wilcoxonCI(pred,target, alpha)
%	
%	See also WILCOXONEXACT, WILCOXONVARIANCE

%	References:
%       Sen 1960
%       DeLong et al, Biometrics, September 1988
%	

%	Copyright 2011 Alistair Johnson

%   $LastChangedBy: alistair $
%   $LastChangedDate: 2014-07-30 10:31:51 +0100 (Wed, 30 Jul 2014) $
%   $Revision: 1182 $
%   Originally written on GLNXA64 by Alistair Johnson, 18-Oct-2011 12:00:52
%   Contact:alistairewj@gmail.com


% First parse the inputs into matrices X and Y
%   X: predictions for target=1
%   Y: predictions for target=0
% Each prediction will be stored in a column of X and Y


if nargin<3
    alpha=0.05;
end
if nargin<2
    error('Incorrect number of input arguments.')
end

% Parse first two inputs to get sizes (m,n) and vectors (X,Y)
K=1; % Only one prediction/target pair
idx=target==1;
X=pred(idx); m=length(X);
Y=pred(~idx); n=length(Y);

% % % -------------- AUC, V10, V01 -------------- % % %
% Using matrices X and Y, calculate estimated Wilcoxon statistic (auc)
% Also Calculate the mxK and nxK V10 and V01 matrices

auc=zeros(1,K);
V10=zeros(m,K); 
V01=zeros(n,K);

for r=1:K % For each X/Y column pair
    % compare 0s to 1s
    for i=1:m
        phi1=sum(gt(X(i,r),Y(:,r))); % Xi>Y
        phi2=sum(eq(X(i,r),Y(:,r))); % Xi=Y
        V10(i,r)=(phi1+phi2*0.5)/n;
        auc(r)=auc(r)+phi1+phi2*0.5;
    end
    auc(r)=auc(r)/(n*m);
    for j=1:n
        phi1=gt(X(:,r),Y(j,r)); % X>Yj
        phi2=eq(X(:,r),Y(j,r)); % X=Yj
        V01(j,r)=(sum(phi1)+sum(phi2)*0.5)/m;
    end
end

%  Calculate S01 and S10, covariance matrices of V01 and V10
S01 = ((V01'*V01)-n*(auc'*auc))/(n-1);
S10 = ((V10'*V10)-m*(auc'*auc))/(m-1);
% S01a = ((V01-auc)'*(V01-auc))/(n-1); % Alternative equivalent formulation
% S10a = ((V10-auc)'*(V10-auc))/(m-1); % Alternative equivalent formulation


% Combine for S, covariance matrix of auc
S = (1/m)*S10 + (1/n)*S01;
auc_ci = norminv([alpha/2,1-(alpha/2)],auc,sqrt(S));
end

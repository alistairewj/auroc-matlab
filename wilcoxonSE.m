function [ se ] = wilcoxonSE(W,Q1,Q2,N0,N1) 
%WILCOXONSE	Calculate the standard error of a Wilcoxon statistic
%	[ se ] = wilcoxonSE(W,Q1,Q2,N0,N1)  calculates the standard error of 
%	the Wilcoxon statistic using probability values Q1 and Q2, the number
%	of negative outcomes N0, and the number of positive outcomes N1. 
%
%	Inputs:
%		W - Probability( PRED|y=1 > PRED|y=0 ) ... AKA Wilcoxon statistic
%		Q1 - Probability( Two PRED|y=1 > One PRED|y=0 )
%       Q2 - Probability( One PRED|y=1 > Two PRED|y=0 )
%           ... where y is the target variable, either 0 or 1, and PRED is
%           the prediction made (usually the probability that y=1).
%		N0 - Number of observations with y=0
%		N1 - Number of observations with y=1
%
%	Outputs:
%		se - Standard error associated with W
%
%	Example
%       pred=rand(100,1); 
%       target=round(pred);
%       target(1:10:end)=1-target(1:10:end);
%       W = wilcoxonEXACT(pred,target);
%       [Q1,Q2] = wilcoxonQ(pred,target);
%       N0=sum(target==0);
%       N1=sum(target==1);
%		[ se ] = wilcoxonSE(W,Q1,Q2,N0,N1) 
%	
%	See also FCN1

%	References:
%       Hanley1982
%	
%	

%	Copyright 2011 Alistair Johnson

%   $LastChangedBy: alistair $
%   $LastChangedDate: 2014-07-30 10:31:51 +0100 (Wed, 30 Jul 2014) $
%   $Revision: 1182 $
%   Originally written on GLNXA64 by Alistair Johnson, 17-Oct-2011 12:11:42
%   Contact:alistairewj@gmail.com

W2=W^2;
se=sqrt((W-W2+(N1-1)*(Q1-W2)+(N0-1)*(Q2-W2))/(N0*N1));

end
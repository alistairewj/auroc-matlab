function [ Q1, Q2, W ] = wilcoxonQ(i1,i2)
%WILCOXONQ	Calculate the Q1 and Q2 probabilities associated with the Wilcoxon statistic.
%
%	[ Q1, Q2 ] = wilcoxonQ( W ) calculates the Q1 and Q2 probabilities
%	given the Wilcoxon statistic (AUROC)
%	
%   [ Q1, Q2 ] = wilcoxonQ( pred, target) calculates the Q1 and Q2
%   probabilities from the Wilcoxon statistic calculated from the inputs
%   pred (predicted probabilities) and target (target binary values)
%
%   [ Q1, Q2, W ] = wilcoxonQ(...) also returns the Wilcoxon statistic.
%
%       Note that this program uses formulas provided in Hanley1982, which
%       assume that the distribution of the area under the receiver
%       operator characteristic (AUROC), or equivalently the Wilcoxon
%       statistic, is a negative exponential. This is the most conservative
%       distributional assumption when compared to Gaussian or gamma.
%
%	Inputs:
%       pred   - the prediction made (the probability that target=1).
%       target - the target variable, 0 or 1, predicted by pred
%       W      - AUROC // Wilcoxon statistic
%
%	Outputs:
%		Q1 - Probability( Two PRED|target=1 > One PRED|target=0 )
%       Q2 - Probability( One PRED|target=1 > Two PRED|target=0 )
%
%	Example
%       pred=rand(100,1); 
%       target=round(pred);
%       target(1:10:end)=1-target(1:10:end);
%       W = wilcoxonEXACT(pred,target);
%       [Q1,Q2] = wilcoxonQ(W)
%       [Q1,Q2] = wilcoxonQ(pred,target)
%	
%	See also WILCOXONQ WILCOXONEXACT WILCOXONSE

%	References:
%       Hanley1982

%   $LastChangedBy: alistair $
%   $LastChangedDate: 2014-07-30 10:31:51 +0100 (Wed, 30 Jul 2014) $
%   $Revision: 1182 $
%   Originally written on GLNXA64 by Alistair Johnson, 17-Oct-2011 12:20:44
%   Contact:alistairewj@gmail.com

if nargin==2
    pred=i1;
    target=i2;
    W=wilcoxonEXACT(pred,target); % Get AUROC
elseif nargin==1
    if length(i1)==1 % Scalar
        % Assume 1st input is AUC
        W=i1;
    else
        error('If only one argument is provided, it must be a scalar (AUROC).\n');
    end
else
    error('Incorrect input format. See help.\n');
end

Q1=W/(2-W);
Q2=2*W*W/(1+W);

end
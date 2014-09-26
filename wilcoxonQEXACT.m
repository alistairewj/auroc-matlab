function [ Q1, Q2, W ] = wilcoxonQEXACT(pred,target)
%WILCOXONQ	Calculate the Q1 and Q2 probabilities associated with the
%   Wilcoxon statistic. It makes no distributional assumptions.
%
%   [ Q1, Q2 ] = wilcoxonQEXACT( pred, target) calculates the Q1 and Q2
%   probabilities from the Wilcoxon statistic calculated from the inputs
%   pred (predicted probabilities) and target (target binary values)
%
%   [ Q1, Q2, W ] = wilcoxonQEXACT(...) also returns the Wilcoxon statistic
%   
%	Inputs:
%       pred   - the prediction made (the probability that target=1).
%       target - the target variable, 0 or 1, predicted by pred
%
%	Outputs:
%		Q1  - Probability( Two PRED|target=1 > One PRED|target=0 )
%       Q2  - Probability( One PRED|target=1 > Two PRED|target=0 )
%       W   - The Wilcoxon statistic
%
%	Example
%       pred=rand(100,1);
%       target=round(pred);
%       target(1:10:end)=1-target(1:10:end);
%       [Q1,Q2,W] = wilcoxonQEXACT(pred,target)
%
%	See also WILCOXONQ WILCOXONEXACT WILCOXONSE

%	References:
%       Hanley1982

%	Copyright 2011 Alistair Johnson

%   $LastChangedBy: alistair $
%   $LastChangedDate: 2014-07-30 10:31:51 +0100 (Wed, 30 Jul 2014) $
%   $Revision: 1182 $
%   Originally written on GLNXA64 by Alistair Johnson, 17-Oct-2011 12:20:44
%   Contact:alistairewj@gmail.com

if nargin~=2
    error('This function requires two inputs (predictions and targets).\n');
end

W=wilcoxonEXACT(pred,target); % Get AUROC

idx=target==0;
Y=sort(pred(idx),'ascend'); % Predictions for negative outcomes
X=sort(pred(~idx),'ascend'); % Predictions for positive outcomes
N0=length(Y);
N1=length(X);

% Q1: Probability two randomly chosen positives are ranked higher than a
% randomly chosen negative
% Q2: Probability one randomly chosen positive is ranked higher than two
% randomly chosen negatives
Q1=0;
Q2=0;
for m=1:N1-1 % For all positive probabilities
    % Find first location of equivalent negative probability
    % This is a fast implementation using histc, taking advantage of the
    % fact that Y contains monotonically non-decreasing values
    [tmpIdx,tmpIdx]=histc(X(m),Y);
    if tmpIdx==0 % X(m) not found within Y
        if X(m)<Y(1)
            tmpIdx=1; % X(m) is below lowest Y value
        elseif X(m)>Y(end)
            tmpIdx=N0; % X(m) is above highest Y value
        else
            % This statement should never execute
            fprintf('Warning: X(m) not found within Y, and not greater or less than the max/min values.\n');
        end
    else
        if abs(X(m)-Y(tmpIdx))<abs(X(m)-Y(tmpIdx+1))
            % Do not change tmpIdx
        else
            tmpIdx=tmpIdx+1; % Increment tmpIdx
        end
    end
    % Now tmpIdx has location of closest probability, tmpVal, to current
    % positive probability X(m)
    
    % For Q1, all values below X(m) are also below X(m+1):X(end)
    % Therefore, there are tmpIdx Y values less than X(m), which has N1-m 
    % values it could be combined with...X(end)
    Q1=Q1+tmpIdx*(N1-m-1);
    
    % For Q2, all n values below Y(tmpIdx) are below X(m)
    % How many combinations of Y(1:tmpIdx) are there below X(m)?
    % Can choose tmpIdx*(tmpIdx-1) + (tmpIdx-1)*(tmpIdx-2) + ... + 2*1
    % combinations of Y values
    % This is equivalent to nchoosek(tmpIdx,2)
    Q2=Q2+nchoosek(tmpIdx,2);
end

% Hypothetically, the maximum number of pairs for Q1 is: N0*nchoosek(N1,2)
% Similarly, the maximum number of pairs for Q2 is: N1*nchoosek(N0,2)
Q1=Q1/(N0*nchoosek(N1,2));
Q2=Q2/(N1*nchoosek(N0,2));

end
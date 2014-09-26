function [ thetaP, theta2 ] = wilcoxonConfidence(L, S, theta, alpha )
%WILCOXONCONFIDENCE	Calculates the covariance for R sets of predictions 
%   and outcomes
%
%	[ thetaP ] = wilcoxonConfidence(L, S, theta ) calculates, given an
%	appropriate matrix of contrasts L, whether the desired AUROC is
%	signficantly better than any of the other methods
%
%	[ thetaP, thetaCI ] = wilcoxonConfidence(L, S, theta ) also returns
%	confidence intervals for theta
%
%	[ thetaP, X2stat ] = wilcoxonConfidence(L, S, theta ) also returns the
%	chi2 statistic
%
%	[ ... ] = wilcoxonConfidence(L, S, theta, alpha ) allows the user
%	to specify the critical value alpha (default 0.05)
%
%	Inputs: 
%       L       - The contrast matrix (rows must sum to 0)
%		S       - Covariance matrix of the AUROCs
%       theta   - Column vector containing AUROC/Wilcoxon statistics
%       alpha   - Critical value used in confidence intervals (default 0.05)
%
%	Outputs:
%		thetaP  - P-value indicating significance
%       thetaCI - Confidence intervals
%		X2stat  - Chi2 statistic for significance testing
%
%	Example:
%   If one wants to determine if a single model is signficantly better than
%   the others, use a contrast row vector with 1 for the desired model and
%   -1/(R-1) for the other models, where R is the total number of
%   models.
%   
%   For three classifiers:
%   [ S, S10, S01, V10, V01, theta ] = ...
%           wilcoxonCovariance(pred1, target1, pred2, target2, pred3, target3) 
%   
%   To test if there is improvement for model 2 over the others, use
%   contrast vector:
%   
%   L=[-1/2 1 -1/2];
%	[ thetaP, thetaCI ] = ...
%           wilcoxonConfidence(L, S, theta, 0.05 ) 
%   
%	thetaP - The p-value of significance for the 2nd classifier
%	thetaCI - Confidence intervals on the 2nd AUROC
%
%	The difference can then be evaluated by the confidence intervals
%
%   To test if model 2 is better than at least one of the other models:
%   
%   L=[-1,1,0;0,1,-1];
%	[ thetaP, X2stat ] = ...
%           wilcoxonConfidence(L, S, theta, 0.05 ) 
%   
%	See also WILCOXONEXACT, WILCOXONVARIANCE

%	References:
%       DeLong et al, Biometrics, September 1988
%	

%	Copyright 2011 Alistair Johnson

%   $LastChangedBy: alistair $
%   $LastChangedDate: 2014-08-12 18:45:11 +0100 (Tue, 12 Aug 2014) $
%   $Revision: 1186 $
%   Originally written on GLNXA64 by Alistair Johnson, 18-Oct-2011 12:00:52
%   Contact:alistairewj@gmail.com

if nargin<3
    error('Not enough input arguments.\n');
end
if nargin<4
    alpha=0.05;
end

if alpha < 0 || alpha > 1
    error('alpha (critical ratio) must be between 0 and 1 (default 0.05).');
elseif alpha > 0.5
    alpha = 1 - alpha;
end

L_sz=size(L);
S_sz=size(S);
theta_sz=size(theta);

% Ensure thetas are a column vector
if theta_sz(1)==1
    theta=theta'; % Transpose to column
elseif theta_sz(2)==1
    % Do nothing, theta is a column vector
else
    % Theta is a matrix - error
    error('Input theta should be a column vector.');
end

if S_sz(2)~=L_sz(2) % Contrast column # is not equal to number of classifiers
    error('Dimension mismatch between contrast L and covariance S.');
end
LSL=L*S*L';
if L_sz(1)==1 % One row
    % Compute using the normal distribution
    mu=L*theta;
    sigma=sqrt(LSL);
    % theta1 = normpdf(0,mu,sigma);
    thetaP = normcdf(0,mu,sigma); 
    % 2-sided test, double the tails -> double the p-value
    if mu<0
        thetaP=2*(1-thetaP);
    else
        thetaP=2*thetaP;
    end
    theta2=norminv([alpha/2,1-alpha/2],theta(L==1),sigma);
else
    % Calculate chi2 stat with DOF=rank(L*S*L');
    w_chi2=theta'*L'*(inv(LSL))*L*theta;
    w_df=rank(LSL);
    thetaP=1-chi2cdf(w_chi2,w_df);
    theta2=w_chi2;
end

end
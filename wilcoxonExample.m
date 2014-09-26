%WILCOXONEXAMPLE	
%	Example running the Wilcoxon code with real predictions from a dataset
%	of critically ill intensive care unit patients.

%	References:
%	The predictions used here were calculated on the Physionet/CinC 2012
%	challenge database, set-a. See: www.physionet.org/challenge/2012
%	This database is freely available to the public.
%	
%	Physionet is described here:
%	Golderberger, 2000.
%
%	Good references on Wilcoxon/AUROC analysis include Delong, 1988 and
%	Song, 1997.

%	Copyright 2014 Alistair Johnson

%	$LastChangedBy$
%	$LastChangedDate$
%	$Revision$
%	Originally written on GLNXA64 by Alistair Johnson, 30-Jul-2014 11:37:23
%	Contact: alistairewj@gmail.com


%=== First load the data
fp = fopen('wilcoxonexample.csv','r');
header = fgetl(fp); header = regexp(header,',','split');
pred = textscan(fp,'%f%f%f%f%f','delimiter',',');
fclose(fp);

target = pred{1};
header = header(2:end);
pred = horzcat(pred{2:end});
P = size(pred,2);

%=== Calculate the AUROC using the formula for the Wilcoxon statistic
%*** NOTE: This version does NOT count ties!
[ W ] = wilcoxon( pred, target );

%=== Calculate the AUROC using the formula for the Wilcoxon statistic which
%factors in ties
[ Wexact ] = wilcoxonEXACT( pred, target );

% Also note that the simple method is much faster due to vectorization.
fprintf('AUROC for various models:\n');
for p=1:P
    fprintf('\t%s',header{p});
end
fprintf('\n');


% Print out the approximate AUROC - this excludes ties
fprintf('Approx');
for p=1:P
    fprintf('\t%0.3f',W(p));
end
fprintf('\n');

% Print out the exact AUROC - this includes ties
fprintf('Exact');
for p=1:P
    fprintf('\t%0.3f',Wexact(p));
end
fprintf('\n');



%% Calculate covariance
[ S, S10, S01, V10, V01, theta ] = wilcoxonCovariance(pred,target);

% The main diagonal of S contains the variance of the AUROC
fprintf('\tAUROC\tVar\t95%% CI\n');
for p=1:P
    % Calculate confidence interval as 2*(Standard Error)
    tmp_CI = 2*sqrt(S(p,p));
    tmp_CI = [theta(p) - tmp_CI, theta(p) + tmp_CI];
    fprintf('%s\t%0.3f\t%0.3f\t[%0.3f,%0.3f]\n',...
        header{p},theta(p),S(p,p),tmp_CI(1),tmp_CI(2));
end
fprintf('\n');
% Note that these estimates do NOT include the correlation between predictions
% Thus even though they overlap, some predictors may be significantly better

%% Is the SVM better than the NN?
% Define a contrast which compares the SVM (2nd model) against NN (4th model)
L = [0,1,0,-1]; % must always sum to 0
alpha = 0.05; % significance level
[ thetaP, thetaCI ] = wilcoxonConfidence(L, S, theta, alpha );

fprintf('Is SVM (2nd model) significantly different than NN (4th model)?\n');
fprintf('SVM (%0.3f) vs NN (%0.3f) - p = %0.4f\n',theta(2),theta(4),thetaP);
fprintf('95%% CI centered on SVM: [%0.3f,%0.3f]\n',thetaCI(1),thetaCI(2));

if thetaP<alpha
    fprintf('Yes, it is significantly different, p < %0.2f\n',alpha);
else
    fprintf('No, it is not significantly different, p > %0.2f\n',alpha);
end
fprintf('\n');

%% Is the SVM better than the RLR?
% Define a contrast which compares the SVM (2nd model) against RLR (1st model)
L = [-1,1,0,0]; % must always sum to 0
alpha = 0.05; % significance level
[ thetaP, thetaCI ] = wilcoxonConfidence(L, S, theta, alpha );

fprintf('Is SVM (2nd model) significantly different than RLR (1st model)?\n');
fprintf('SVM (%0.3f) vs RLR (%0.3f) - p = %0.4f\n',theta(2),theta(1),thetaP);
fprintf('95%% CI centered on SVM: [%0.3f,%0.3f]\n',thetaCI(1),thetaCI(2));

if thetaP<alpha
    fprintf('Yes, it is significantly different, p < %0.2f\n',alpha);
else
    fprintf('No, it is not significantly different, p > %0.2f\n',alpha);
end
fprintf('\n');

%% Is X better than Y?
% The following code allows you to play around with different comparisons
cmp1 = 1;
cmp2 = 4;
alpha = 0.05; % significance level


L = zeros(1,4); % must always sum to 0
L(cmp1) = 1;
L(cmp2) = -1;
[ thetaP, thetaCI ] = wilcoxonConfidence(L, S, theta, alpha );

fprintf('Is %s significantly different than %s?\n', header{cmp1}, header{cmp2});
fprintf('%s (%0.3f) vs %s (%0.3f) - p = %0.4f\n',header{cmp1},theta(2),header{cmp2},theta(1),thetaP);
fprintf('95%% CI centered on %s: [%0.3f,%0.3f]\n',header{cmp1},thetaCI(1),thetaCI(2));

if thetaP<alpha
    fprintf('Yes, it is significantly different, p < %0.2f\n',alpha);
else
    fprintf('No, it is not significantly different, p > %0.2f\n',alpha);
end
fprintf('\n');


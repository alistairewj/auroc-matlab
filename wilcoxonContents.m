% wilcoxonContents
% Displays info regarding various wilcoxon functions
% Copyright 2011 Alistair Johnson

fprintf('[ W ] = wilcoxon( pred, target )\n');
fprintf('Calculates AUROC using vectorized sorting algorithm\n');

fprintf('\n'); % Spacing

fprintf('[ W ] = wilcoxonEXACT( pred, target )\n');
fprintf('Calculates AUROC using case comparison algorithm\n');
fprintf('(adds 0.5 in cases of a tie). \n');
fprintf('These are different functions because wilcoxonEXACT is simpler, \n');
fprintf('but takes an order of magnitude more time to compute.\n');

fprintf('\n'); % Spacing

fprintf('[ W ] = wilcoxonEXACTXY( X, Y )\n');
fprintf('Identical to wilcoxonEXACT, except uses vectors X and Y,\n');
fprintf('containing only positive and negative predictions, respectively.\n');

fprintf('\n'); % Spacing

fprintf('[ S, S10, S01, V10, V01, theta ] = ...\n');
fprintf('wilcoxonCovariance(pred1,target1,pred2,target2,...)\n');
fprintf('Calculates covariance matrices for multiple pred/target pairs.\n');

fprintf('\n'); % Spacing

fprintf('[ Q1, Q2 ] = wilcoxonQ( W )\n');
fprintf('Calculates the Q1 and Q2 probabilities given the Wilcoxon statistic (AUROC)\n');

fprintf('\n'); % Spacing

fprintf('[ Q1, Q2 ] = wilcoxonQEXACT( pred, target ) \n');
fprintf('Calculates the Q1 and Q2 probabilities from the \n');
fprintf('Wilcoxon statistic calculated from the inputs \n');
fprintf('pred (predicted probabilities) and target (target binary values)\n');

fprintf('\n'); % Spacing

fprintf('[ se ] = wilcoxonSE(W,Q1,Q2,N0,N1)\n');
fprintf('Calculates the standard error of the Wilcoxon statistic \n');
fprintf('using probability values Q1 and Q2, the number of negative outcomes N0, \n');
fprintf('and the number of positive outcomes N1. \n');

fprintf('\n'); % Spacing

fprintf('[ thetaCI ] = wilcoxonConfidence(L, S, theta ) \n');
fprintf('Calculates confidence intervals using the contrast row vector L \n');

fprintf('\n'); % Spacing

fprintf('[ theta, thetaCI, S, S10, S01, V10, V01 ] = wilcoxonCI(pred,target,alpha)\n');
fprintf('Calculates the variance and confidence interval for a single set of\n');
fprintf('predictions and outcomes\n');


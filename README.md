auroc-matlab
============

This repository facilitates the calculation of the Area Under the Receiver Operator Characteristic curve, or AUROC. This is a commonly used metric for assessing the discrimination of a set of predictions for a binary target (i.e. for predicting a 0 or a 1). While traditionally used with probabilistic predictions (i.e. predictions lying within the 0-1 interval), the AUROC is a rank based metric and does not require predictions to lie within the 0-1 interval. The AUROC is related to the Mann-Whitney U statistic (up to a normalising factor).

This repository also implements the non-parametric statistical significance test as published by DeLong and DeLong (1988). This test allows for the comparison of two AUROCs calculated on the same dataset. There will be a strong correlation for AUROCs calculated on the same dataset, and the statistical test used must account for this. The technique by DeLong and DeLong is one method - another is an empirical adjustment proposed by Hanley and McNeil (1982, 1983) though these are not implemented here.

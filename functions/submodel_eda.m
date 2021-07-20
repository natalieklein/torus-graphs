function [marginal,phase_diffs,phase_sums]=submodel_eda(X)
% This function uses the CircStat2012a toolbox, 
% freely available at: 
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

% input
% X, d by N, matrix of angle vector of observations over repeated trials

% output 
% individual p-values testing the concentrations of marginals, phase
% differences, and phase sums. 

d=size(X,1);
%% testing for uniformity in marginals
% null hypothesis is that the concentration is circular uniform
pval_samp_marg=nan(1,d);
for ind=1:d
    pval_samp_marg(ind)= circ_rtest(X(ind,:));
end
pval_samp_marg(pval_samp_marg==0)=pval_samp_marg(pval_samp_marg==0) + 1e-100; % in case some pvals come out as zero, to prevent log(0)=-inf

% combining the p-vals in a fisher test
Comb_stat_samp_marg = -2*sum(log(pval_samp_marg));
dofs=2*length(pval_samp_marg);
group_pval_samp_marg = 1 - chi2cdf(Comb_stat_samp_marg,dofs); 

marginal.pvals=pval_samp_marg;
marginal.group_pval=group_pval_samp_marg;

%% testing for uniformity in phase sums and differences
num_pairs=d*(d-1)/2;
pval_samp_sum=nan(1,num_pairs);
pval_samp_diff=nan(1,num_pairs);
inc=1;
for ind_i=1:d
    for ind_j=(ind_i+1):d
        pval_samp_sum(inc)= circ_rtest(wrapToPi(X(ind_i,:)+X(ind_j,:)));
        pval_samp_diff(inc)= circ_rtest(wrapToPi(X(ind_i,:)-X(ind_j,:)));
        inc=inc+1;
    end
end


pval_samp_diff(pval_samp_diff==0)=pval_samp_diff(pval_samp_diff==0) + 1e-100; % in case some pvals come out as zero, to prevent log(0)=-inf
pval_samp_sum(pval_samp_sum==0)=pval_samp_sum(pval_samp_sum==0) + 1e-100; % in case some pvals come out as zero, to prevent log(0)=-inf

% combining the p-vals in a fisher test
Comb_stat_samp_sum = -2*sum(log(pval_samp_sum));
Comb_stat_samp_diff = -2*sum(log(pval_samp_diff));
dofs=2*num_pairs;

group_pval_samp_sum = 1 - chi2cdf(Comb_stat_samp_sum,dofs); 
group_pval_samp_diff = 1 - chi2cdf(Comb_stat_samp_diff,dofs);


phase_sums.pvals=pval_samp_sum;
phase_sums.group_pval=group_pval_samp_sum;

phase_diffs.pvals=pval_samp_diff;
phase_diffs.group_pval=group_pval_samp_diff;
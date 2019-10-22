function [FPR,FNR]=false_rates(sim_true_edges,p_values,alpha_level)
% sim_true_edges is a logical vector (length is num possible edges), where
% 1 means edge, and 0 means no edge
% pvals is num possible edges by num simulations

% output
% scalar for FPR and FNR
num_edges=length(sim_true_edges);
%% We are interested in the FPR and FNR
P=sum(sim_true_edges);%true positives
N=num_edges-P;%true negatives
pval_decision=p_values < alpha_level;

inds_rho1=sim_true_edges;
inds_rho0=~sim_true_edges;

% Considering each simulation separately
TPR_sims=sum(pval_decision(inds_rho1,:),1)/P;
FPR_sims=sum(pval_decision(inds_rho0,:),1)/N;

TPR=mean(TPR_sims); % true positive (hits),
FPR=mean(FPR_sims);% false positive (false alarm),

FNR=1-TPR;


%Output: Investigation of the ability of torus graphs to recover the true structure 
% as a function of true edge density, sample size, and data dimension
%Related figures: Fig. 5 (main text) and Fig. S4 (supplements)

% it takes a long time to run. 
clear all;%close all;clc
addpath(genpath('functions'))
rng(10) % for reproducibility
k=0; % 0 marginal concentration
n_list=[200 400 840 1600 3200];%Sample size
all_prop_edges_on=[0.25 0.5]; % edge density
d_all=[4,8,16,24,32]; % number of dimensions
num_edges_all=d_all.*(d_all-1)/2;

num_sim=30; % number of simulations to run in each case

num_props=length(all_prop_edges_on);
num_nlist=length(n_list);
num_d=length(d_all);
all_pvals=cell(1,num_d);

optO.burnin=200;optO.thin=50;
%%
for ni=1:num_nlist % Loop through sample size
    optO.nsamp=n_list(ni);
    for di=1:num_d % loop through number of dimensions
        d=d_all(di);
        num_edges=num_edges_all(di);
        for pri=1:num_props % loop through edge density
            disp([ni/num_nlist,di/num_d,pri/num_props])
            prop_edges_on=all_prop_edges_on(pri)*num_edges;
            if mod(prop_edges_on,1) == 0 % checking that we have an integer number for edges
                temp_edges=zeros(num_edges,1);
                temp_edges(1:prop_edges_on)=1;
                pval_allr=nan(num_edges,num_sim);
                parfor r=1:num_sim % loop through simulations for each case
                    % generating samples 
                    phi_true=[  k*cos(zeros(d,1));...
                                k*sin(zeros(d,1));...
                                temp_edges;...
                                zeros(num_edges,1);...
                                zeros(num_edges,1);...
                                zeros(num_edges,1)];
                    Xsim_samples=sampleGibbs_tutorial(d, phi_true, optO);
                    
                    % fitting full torus graphs, though the underlying
                    % model is only phase differences with zero marginals
                    [~,edges]=torus_graphs(Xsim_samples);
                    % fitting torus graphs with underlying model
%                    [~,edges]=torus_graphs(Xsim_samples,[],[],[false true false]);
                    
                    pval_allr(:,r)=edges.p_vals;
                end
                all_pvals{di}(:,:,pri,ni)=pval_allr;
            end
        end
    end
end
%% varying the p-value thersholds to get edge decision for ROC curves
% note that thresholding the p-value is equivalent to varying the threshold
% on the edgewise chi-squared statistic.
alpha_thresh_list=linspace(1e-200,1-1e-200,20e3);
num_thresh=length(alpha_thresh_list);
AUC_all=nan(num_d,num_nlist,num_props);
TPR_pval=nan(num_thresh,num_d,num_nlist,num_props);
FPR_pval=nan(num_thresh,num_d,num_nlist,num_props);
PRE_pval=nan(num_thresh,num_d,num_nlist,num_props);
for ni=1:num_nlist
    for di=1:num_d
        num_edges=num_edges_all(di);
        for pri=1:num_props
            prop_edges_on=all_prop_edges_on(pri)*num_edges;
            P=prop_edges_on;%true positives
            N=num_edges-P;%true negatives
            if mod(prop_edges_on,1) == 0 % checking that we have whole edges
                pval_dec=nan(num_thresh,num_edges,num_sim);
                for i=1:num_thresh
                    alpha_thresh=alpha_thresh_list(i);
                    pval_dec(i,:,:)=all_pvals{di}(:,:,pri,ni) <= alpha_thresh;
                end
                temp_edges=false(num_edges,1);
                temp_edges(1:prop_edges_on)=true;
                inds_rho1=temp_edges;
                inds_rho0=~temp_edges;
                % computing ROC curve
                TP=sum(sum(pval_dec(:,inds_rho1,:),3),2);
                FP=sum(sum(pval_dec(:,inds_rho0,:),3),2);
                
                TPR_pval(:,di,ni,pri)=TP/(P*num_sim); % true positive (hits), 101 to 200 sims are truly one
                FPR_pval(:,di,ni,pri)=FP/(N*num_sim);% false positive (false alarm), 1 to 100 sims are truly zero
                
                PRE_pval(:,di,ni,pri)=TP./(TP+FP);% TP/(TP+FP)
                AUC_all(di,ni,pri)=trapz(FPR_pval(:,di,ni,pri),TPR_pval(:,di,ni,pri));
                
            end
        end
    end
    
    figure('color','w');
    for pri=1:num_props
        subplot(2,3,pri);hold on
        plot(FPR_pval(:,:,ni,pri),TPR_pval(:,:,ni,pri),'linewidth',2)
        xlabel('FPR'); ylabel('TPR')
        legend({'4','8','16', '24', '32'})
        title(sprintf('n=%i',n_list(ni)))
        set(gca,'fontsize',14)
        
        subplot(2,3,pri+3);hold on
        plot(TPR_pval(:,:,ni,pri),PRE_pval(:,:,ni,pri),'linewidth',2)
        xlabel('Recall (TPR)');ylabel('Precision (1-FDR)')
        legend({'4','8','16', '24', '32'})
        title(sprintf('n=%i',n_list(ni)))
        set(gca,'fontsize',14)
    end
    
end
%%
figure('color','w');
for pri=1:num_props
    subplot(1,3,pri);hold on
    AUC_all(di,ni,pri)
    for di=1:num_d,plot(n_list,AUC_all(di,:,pri),'o-','linewidth',2),end
    xlabel('Sample size'); ylabel('ROC AUC')
    legend({'4','8','16', '24', '32'})
    set(gca,'fontsize',14)
    set(gca,'xscale','log')
    set(gca,'xtick',n_list,'xticklabel',n_list,'XMinorGrid','off','XMinorTick','off')
    ylim([0.4 1])
    grid on
end


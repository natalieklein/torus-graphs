% Start here; this script shows you how to use the torus_graphs function 
% Output: results for Torus graph analysis of coupling in 24-dimensional LFP data with four distinct regions
% Related figures: Fig. 8 (main text) and S8 (supplements)
load('data/anglebeta.mat')
addpath(genpath('functions'))
X=anglebeta; % phase angle data with 24 dimensions and 840 trials
% we have labels for the subregions:
sub_reg_inds.ca3=strcmp(subregions,'CA3');
sub_reg_inds.dg=strcmp(subregions,'DG');
sub_reg_inds.sub=strcmp(subregions,'Sub');
sub_reg_inds.pfc=strcmp(subregions,'PFCv');
%% Quick implementation, fit TGs to obtain a 24dim network, Fig 8B,C
TG=torus_graphs(X);
figure;
subplot(121);network_graph=plot(TG); % to plot the network graph
% customizing the location of the nodes in the network plot
custom_node_placement(network_graph,sub_reg_inds)
subplot(122);imagesc(full(TG.adjacency)'); 
% In, the text we ignore the presence of the edge between DG and PFC
% because Fig8A tests did not indicate evidence for cross-region edges
% between DG and PFC. 
%% Choosing a sub-model of TG
% we can test the concentrations of marginals, and pairwise phase
% differences and sums. 
%The null hypothesis is that there is no circular concentration
[marginal,phase_diffs,phase_sums]=submodel_eda(X)
% for this dataset we reject the null hypothesis of null concentrations in
% phase sums. So we decide to fit the phase difference model
sel_mode=[true true false]; % entries in this vector correspond to [marginal, phase differences, phase sums]
% This submodel only considers marginal and phase difference concentrations
TG_sub=torus_graphs(X,[],[],sel_mode);
figure;
subplot(121);network_graph_submodel=plot(TG_sub); % to plot the network graph
custom_node_placement(network_graph_submodel,sub_reg_inds)
subplot(122);imagesc(full(TG_sub.adjacency)'); 

% note that for this dataset the resulting graph is almost identical to
% that of the full model [true true true] (default)
%% Grouping electrodes by region and choosing alpha level, Fig. 8A
alpha_level=0.001/6;%Bonferroni, 6 possible edges. Default is just 0.05
dj=[5,8,3,8]; % grouping information, numer of electrodes in each region
sel_mode=[true true false];
[TG_regions,edges,phi_hat,inference]=torus_graphs(X,alpha_level,dj,sel_mode);
figure;
network_graph_regions=plot(TG_regions,'NodeLabel',{'CA3','DG','Sub','PFC'}); % to plot the network graph
network_graph_regions.XData(1)=21;network_graph_regions.YData(1)=5; %ca3
network_graph_regions.XData(2)=2;network_graph_regions.YData(2)=6;%dg
network_graph_regions.XData(3)=-2;network_graph_regions.YData(3)=4; %sub
network_graph_regions.XData(4)=12;network_graph_regions.YData(4)=3;%pfc
axis off
% We see a missing edge between dg and pfc
% the p-values can be inspected here
edges.p_vals
%******************************
% This group analysis uses a stringent alpha level and acts as a overall
% test of cross-region connectivity. Based on the result of no connection
% between dg and pfc, we relax the alpha level for the 24dim network to
% 0.05 and ignore any pairwise-electrode edges between dg and pfc. 
%******************************
%% Sampling from TGs
optO.nsamp=1; % number of samples
optO.burnin=500; % burnin period
optO.thin=100; % thinning
X_sampled=sampleGibbs_tutorial(24,phi_hat,optO);
%% PLV, 24dim network graph, figures S8 
[plv, pvalue, all_edges]=phase_locking_value(X);
 edges_on_plv=all_edges(pvalue<=0.0005,:);
PLV_graph = graph(edges_on_plv(:,1), edges_on_plv(:,2));
figure;
subplot(121);network_graph_PLV=plot(PLV_graph);
custom_node_placement(network_graph_PLV,sub_reg_inds)
subplot(122);imagesc(full(PLV_graph.adjacency)'); 
%% Conditional coupling, case of phase difference model with uniform margins, see Section 4.4. in the paper.
[TG_regions,edges,phi_hat,inference]=torus_graphs(X,0.05,[],[false true false]);
active_edges_inds=edges.p_vals<0.05;
TG_weighted_graph = graph(edges.all_possible(active_edges_inds,1), edges.all_possible(active_edges_inds,2));
scaling_factor=3;
LWidths=scaling_factor*edges.cond_coupling_coeff(active_edges_inds);

figure;
TG_weighted_handle=plot(TG_weighted_graph,'LineWidth',LWidths);
custom_node_placement(TG_weighted_handle,sub_reg_inds)

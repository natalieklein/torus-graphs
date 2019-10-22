%Output: realistic simulated data sets based on positive rotational dependence, 3dim and 5dim 
%Related figures: Fig. 4 (main text), and Fig. S10 (supplements)
rng(20) % for reproducibility
N=840; % sample size
num_sim=200;
%% choose sim type, A (5dim) or B (3dim)
sim_type='A';
% sim_type='B';
if strcmp(sim_type,'A')
    num_edges=5*4/2;%for sim A
    sim_true_edges=logical([1 0 0 0 1 0 0 1 0 1]);
    edge_labels={'1,2','1,3','1,4','1,5','2,3','2,4','2,5','3,4','3,5','4,5'};
elseif strcmp(sim_type,'B')
    num_edges=3*2/2;%for sim B
    sim_true_edges=logical([1 0 1]);
    edge_labels={'1,2','1,3','2,3'};
end
all_pvals_tg=nan(num_edges,num_sim);
all_pvals_plv=nan(num_edges,num_sim);

parfor r=1:num_sim
    Xsim_samples=simAltTG(sim_type,N)';
    
    [~,edges]=torus_graphs(Xsim_samples,[],[],[false true false]);
    all_pvals_tg(:,r)=edges.p_vals;
    
    [~, all_pvals_plv(:,r)]=phase_locking_value(Xsim_samples);
end
%% Results for Figure 4
alpha_level=0.001/num_edges;
[FPR_tg,FNR_tg]=false_rates(sim_true_edges,all_pvals_tg,alpha_level);
[FPR_plv,FNR_plv]=false_rates(sim_true_edges,all_pvals_plv,alpha_level);

figure;
subplot(121);imagesc(all_pvals_tg(:,:)<alpha_level);title(sprintf('TG, FPR=%1.2f, FNR=%1.2f',FPR_tg,FNR_tg))
set(gca,'ytick',1:num_edges,'yticklabel',edge_labels); xlabel('simulations'); ylabel('edges')
subplot(122);imagesc(all_pvals_plv(:,:)<alpha_level);title(sprintf('PLV, FPR=%1.2f, FNR=%1.2f',FPR_plv,FNR_plv))
set(gca,'ytick',1:num_edges,'yticklabel',edge_labels); xlabel('simulations'); ylabel('edges')
%% Results for figure S10
% Visualize simulated and real data to demonstrate the validity of the simulation process.
[Xsim,Xoriginaldata]=simAltTG(sim_type,N);
figure;
subplot(1,2,1);hold on
set(gca,'fontsize',18)
plotmatrix(Xoriginaldata);title('Real data')
subplot(1,2,2);hold on
set(gca,'fontsize',18)
plotmatrix(Xsim);title('Simulated data')

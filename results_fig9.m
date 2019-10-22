%Output: Comparison between phase angles from three LFPs located in PFC and the theoretical torus graph distribution
%Related figures: Fig. 9 (main text)
addpath(genpath('functions'))
load('data/anglebeta.mat')
X=anglebeta;

%% EDA figures, 3vars
Xtrivar=X(21:23,:);
figure(10);subplot(121);plotmatrix(Xtrivar')
%% choosing a submodel
[marginal,phase_diffs,phase_sums]=submodel_eda(Xtrivar)
sel_mode=[true true false];
alpha_level=0.05/3; %Bonferroni
[TG,edges,phi_hat,inference]=torus_graphs(Xtrivar,alpha_level,[],sel_mode);
edges.active
%% evaluate the fitted density
num_steps=50;
step_list=linspace(-pi,pi,num_steps);
[x,y,z]=meshgrid(step_list,step_list,step_list);
f=nan(size(x));
for indx=1:num_steps
    for indy=1:num_steps
        for indz=1:num_steps
            this_x=x(indx,indy,indz);
            this_y=y(indx,indy,indz);
            this_z=z(indx,indy,indz);
            f(indx,indy,indz)=exp(phi_hat(1:6)'*[cos(this_x);cos(this_y);cos(this_z);sin(this_x);sin(this_y);sin(this_z)] +...
                                  phi_hat(7:12)'*[cos(this_x-this_y);cos(this_x-this_z);cos(this_y-this_z);sin(this_x-this_y);sin(this_x-this_z);sin(this_y-this_z)]);
        end
    end
end
norm_const=trapz(step_list,trapz(step_list,trapz(step_list,f,3),2));
f=(1/norm_const)*f;
% let's get pairs
xy_pair=trapz(step_list,f,3);
xz_pair=squeeze(trapz(step_list,f,2));
yz_pair=squeeze(trapz(step_list,f,1));
% let's get marginals
x_marg=trapz(step_list,trapz(step_list,f,3),2);
y_marg=squeeze(trapz(step_list,trapz(step_list,f,3),1));
z_marg=squeeze(trapz(step_list,trapz(step_list,f,2),1));
all_margs{1}=x_marg;all_margs{2}=y_marg;all_margs{3}=z_marg;
%% plot data bivariate histograms vs fitted densities, Fig 9A
figure('color','w');
subplot(331);% using circular histograms for marginals
            h=polarhistogram(Xtrivar(1,:),11,'normalization','pdf'); 
            bin_centers=h.BinEdges(1:end-1)+diff(h.BinEdges)/2;
            bin_values=h.Values;
            cla
            bar(bin_centers,bin_values,'BarWidth',1,'facecolor',[0 .45 .74],'FaceAlpha',.6);hold on
            plot(step_list,x_marg,'linewidth',3,'color',[.85, .33, .1]);axis square     
            axis tight;ylim([0 0.48]);axis off
subplot(332);imagesc(xy_pair);axis xy;axis square;axis off
subplot(333);imagesc(xz_pair);axis xy;axis square;axis off
subplot(334);customBivarHist(Xtrivar([1 2],:),20);axis off
subplot(335);
            h=polarhistogram(Xtrivar(2,:),11,'normalization','pdf'); 
            bin_centers=h.BinEdges(1:end-1)+diff(h.BinEdges)/2;
            bin_values=h.Values;
            cla
            bar(bin_centers,bin_values,'BarWidth',1,'facecolor',[0 .45 .74],'FaceAlpha',.6);hold on
            plot(step_list,y_marg,'linewidth',3,'color',[.85, .33, .1]);axis square
            axis tight;ylim([0 0.48]);axis off
subplot(336);imagesc(yz_pair);axis xy;axis square;axis off
subplot(337);customBivarHist(Xtrivar([1 3],:),20);axis off
subplot(338);customBivarHist(Xtrivar([2 3],:),20);axis off
subplot(339);
            h=polarhistogram(Xtrivar(3,:),11,'normalization','pdf'); 
            bin_centers=h.BinEdges(1:end-1)+diff(h.BinEdges)/2;
            bin_values=h.Values;
            cla
            bar(bin_centers,bin_values,'BarWidth',1,'facecolor',[0 .45 .74],'FaceAlpha',.6);hold on
            plot(step_list,z_marg,'linewidth',3,'color',[.85, .33, .1]);axis square
            axis tight;ylim([0 0.48]);axis off
%% sampling from the fitted model, this portion takes about 2mins
rng(10) % for reproducibility
optO.nsamp=10e3; % number of samples
optO.burnin=500; % burnin period
optO.thin=100; % thinning
tic 
Xsim_trivar=sampleGibbs_tutorial(3,phi_hat,optO);
toc
figure(10);subplot(122);plotmatrix(Xsim_trivar')
%% plot marginal, phase differences and phase sums, Fig 9B
%This portion uses circ_ksdensity to estimate the densities from the simulated data 
% https://www.mathworks.com/matlabcentral/fileexchange/44072-kernel-density-estimation-for-circular-functions?focused=7962020&tab=function
figure('color','w')
d=3;
inc=1;
for i=1:d
    for j=1:d
        subplot(d,d,inc)
        if i==j
            h=polarhistogram(Xtrivar(i,:),11,'normalization','pdf'); 
            bin_centers=h.BinEdges(1:end-1)+diff(h.BinEdges)/2;
            bin_values=h.Values;
            cla
            bar(bin_centers,bin_values,'BarWidth',1,'facecolor',[0 .45 .74],'FaceAlpha',.6);hold on
            plot(step_list,all_margs{i},'linewidth',3,'color',[.85, .33, .1]);
            axis square;axis tight;ylim([0 0.48]);axis off
        else
            if j>i
                sum_example=wrapToPi(Xtrivar(i,:)+Xtrivar(j,:));
                h=polarhistogram(sum_example,11,'normalization','pdf'); 
                bin_centers=h.BinEdges(1:end-1)+diff(h.BinEdges)/2;
                bin_values=h.Values;
                cla
                bar(bin_centers,bin_values,'BarWidth',1,'facecolor',[0.8 .8 .8],'FaceAlpha',.6)
                hold on
                vfPDFSamples=linspace(-pi,pi,20);
                sum_sim_example=wrapToPi(Xsim_trivar(i,:)+Xsim_trivar(j,:));
                vfEstimate=circ_ksdensity(sum_sim_example',vfPDFSamples,[-pi pi],0.4);
                c_lev=1/(2*pi);%mean(vfEstimate);
                plot(vfPDFSamples,c_lev*ones(size(vfPDFSamples)),'linewidth',3,'color',[.85, .33, .1])
                xlim([-pi pi])
                axis square;axis tight;%axis off
                ylim([0 .48])
            else
                diff_example=wrapToPi(Xtrivar(i,:)-Xtrivar(j,:));
                h=polarhistogram(diff_example,11,'normalization','pdf'); 
                bin_centers=h.BinEdges(1:end-1)+diff(h.BinEdges)/2;
                bin_values=h.Values;
                cla
                bar(bin_centers,bin_values,'BarWidth',1,'facecolor',[0.8 .8 .8],'FaceAlpha',.6)
                hold on
                vfPDFSamples=linspace(-pi,pi,200);
                
                diff_sim_example=wrapToPi(Xsim_trivar(i,:)-Xsim_trivar(j,:));
                vfEstimate=circ_ksdensity(diff_sim_example',vfPDFSamples,[-pi pi],0.4);
                plot(vfPDFSamples,vfEstimate,'linewidth',3,'color',[.85, .33, .1])
                xlim([-pi pi])
                axis square;axis tight;%axis off
                ylim([0 0.48])
            end       
        end
%         title([j i])
        axis off
        inc=inc+1;
    end
end
%% ks test for difference of angles
diff_example=wrapToPi(Xtrivar(1,:)-Xtrivar(2,:));
diff_sim_example=wrapToPi(Xsim_trivar(1,:)-Xsim_trivar(2,:));
pValue_ks_diff=[];
[~, pValue_ks_diff(1), ~] = kstest2(diff_example', diff_sim_example');

diff_example=wrapToPi(Xtrivar(1,:)-Xtrivar(3,:));
diff_sim_example=wrapToPi(Xsim_trivar(1,:)-Xsim_trivar(3,:));
[~, pValue_ks_diff(2), ~] = kstest2(diff_example', diff_sim_example');

diff_example=wrapToPi(Xtrivar(2,:)-Xtrivar(3,:));
diff_sim_example=wrapToPi(Xsim_trivar(2,:)-Xsim_trivar(3,:));
[~, pValue_ks_diff(3), ~] = kstest2(diff_example', diff_sim_example');

% Using the Fisher method to combine p-values
pValue_ks_diff
Comb_stat_ks_diff = -2*sum(log(pValue_ks_diff));
dofs=2*length(pValue_ks_diff);% same for both vectors
group_pval_ks_diff = chi2cdf(Comb_stat_ks_diff,dofs,'upper')
%% ks test for marginal concentrations
pValue_ks_marg=[];
[~, pValue_ks_marg(1), ~] = kstest2(Xtrivar(1,:)', Xsim_trivar(1,:)');
[~, pValue_ks_marg(2), ~] = kstest2(Xtrivar(2,:)', Xsim_trivar(2,:)');
[~, pValue_ks_marg(3), ~] = kstest2(Xtrivar(3,:)', Xsim_trivar(3,:)');

pValue_ks_marg
Comb_stat_ks_marg = -2*sum(log(pValue_ks_marg));
dofs=2*length(pValue_ks_marg);% same for both vectors
group_pval_ks_marg = chi2cdf(Comb_stat_ks_marg,dofs,'upper')
%% We cannot reject the null hypothesis that the simulated samples from fitted data 
% come from the same distribution as the original data

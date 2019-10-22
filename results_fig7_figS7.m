%Output: Torus graphs and PLV graphs to infer low-dimensional networks of interest
%Related figures: Fig. 7 (main text), and Fig. S7 (supplements)

%%--------- Fig 7B, S7B: Linear probe in CA3 (5dim)
clear all;load('data/anglebeta.mat')
addpath(genpath('functions'))
ca3_inds=strcmp(subregions,'CA3');
X=anglebeta(ca3_inds,:); 
d=size(X,1);
num_all_possible_edges=d*(d-1)/2;
% PLV estimates and pvalues based on permutation
edges_plv.all_possible=nan(num_all_possible_edges,2);
edges_plv.estimates=nan(num_all_possible_edges,1);
edges_plv.p_vals=nan(num_all_possible_edges,1);
inc=1;
for j=1:d
    for k=(j+1):d
        edges_plv.all_possible(inc,:)=[j,k];
        [edges_plv.estimates(inc),edges_plv.p_vals(inc)]=phase_locking_value(X([j,k],:));
        inc=inc+1;
    end
end
edges_plv.all_possible(edges_plv.p_vals<0.0001,:)
% Torus graphs, pvalues
[~,edges_tg]=torus_graphs(X);
edges_tg.all_possible(edges_tg.p_vals<0.0001,:)
% edges_tg.all_possible(edges_tg.p_vals>=0.0001 & edges_tg.p_vals<0.005,:)
% edges_tg.all_possible(edges_tg.p_vals>=0.0005 & edges_tg.p_vals<0.05,:)

%% --------- Fig 7A, S7A,C: (3dim) combinations: DG, Sub, PFC
clear all;load('data/anglebeta.mat')
dginds = find(ismember(subregions,'DG'));% 6:13
subinds = find(ismember(subregions,'Sub'));% 14:16
pfcvinds = find(ismember(subregions,'PFCv'));%17:24
trivar_combinations = combvec(dginds',subinds',pfcvinds');
% figure; imagesc(trivar_combinations)
% sorting combinations based on DG channels
[~,inds_s]=sort(trivar_combinations(1,:));
trivar_combinations =trivar_combinations(:,inds_s);
% figure; imagesc(trivar_combinations)
num_comb=size(trivar_combinations,2);% 192 trivariate combinations
d=3;
num_all_possible_edges=d*(d-1)/2;
%% PLV estimates and pvalues (takes a couple of minutes)
edges_plv.all_possible=nan(num_all_possible_edges,2);
edges_plv.estimates=nan(num_all_possible_edges,num_comb);
edges_plv.p_vals=nan(num_all_possible_edges,num_comb);
for c=1:num_comb
    sprintf('%i of %i',c,num_comb)
    X=anglebeta(trivar_combinations(:,c),:);
    inc=1;
    for j=1:d
        for k=(j+1):d
            if c==1
                edges_plv.all_possible(inc,:)=[j,k];
            end
            [edges_plv.estimates(inc,c),edges_plv.p_vals(inc,c)]=phase_locking_value(X([j,k],:));
            inc=inc+1;
        end
    end
end
%% Torus graphs, pvalues
edges_tg.p_vals=nan(num_all_possible_edges,num_comb);
for c=1:num_comb
    sprintf('%i of %i',c,num_comb)
    X=anglebeta(trivar_combinations(:,c),:);
    [~,these_edges_tg]=torus_graphs(X);
    edges_tg.p_vals(:,c)=these_edges_tg.p_vals;
end
edges_tg.all_possible=these_edges_tg.all_possible;
%%
figure;
imagesc(edges_plv.p_vals<0.0001)
xlabel('192 trivariate combinations')
title('PLV edges (yellow)')
set(gca,'ytick',1:3,'yticklabel',{'DG-Sub','DG-PFC','Sub-PFC'})

figure;
subplot(311)
imagesc(edges_tg.p_vals<0.0001)
xlabel('192 trivariate combinations')
title('TGs edges (yellow, pval<0.0001)')
set(gca,'ytick',1:3,'yticklabel',{'DG-Sub','DG-PFC','Sub-PFC'})
subplot(312)
imagesc(edges_tg.p_vals>=0.0001 & edges_tg.p_vals<0.005)
xlabel('192 trivariate combinations')
title('TGs edges (yellow, pval in [0.0001,0.005])')
set(gca,'ytick',1:3,'yticklabel',{'DG-Sub','DG-PFC','Sub-PFC'})
subplot(313)
imagesc(edges_tg.p_vals>=0.005 & edges_tg.p_vals<0.05)
xlabel('192 trivariate combinations')
title('TGs edges (yellow, pval in [0.005,0.05])')
set(gca,'ytick',1:3,'yticklabel',{'DG-Sub','DG-PFC','Sub-PFC'})
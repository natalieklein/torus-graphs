function [plv, pvalue, all_edges]=phase_locking_value(Y)
%% simplest case,
% input:
% Y is a 2 by N pair of angles over trials
% output
% plv; phase locking value (scalar) for a pair of angles
% pvalue (scalar) obtained by permutation test
%% Multiple pairs case
% Y is a d by N, where d>2
% output
% plv; phase locking (vector) value for all pairs of angles
% pvalue (vector) obtained using the Rayleigh test
% the Rayleigh test function uses the the CircStat2012a toolbox,
% freely available at:
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

%%
comp_PLV=@(Y) abs(mean(exp(1i*(Y(1,:)-Y(2,:)))));
[d,N]=size(Y);
if d==2
    plv=comp_PLV(Y);
    % obtaining a p-value by permutation
    nperm=20e3;% number of permutations
    PLV_permut=nan(1,nperm);
    parfor p=1:nperm
        random_trials_indices = randsample(N,N);%sample without replacement
        PLV_permut(p)=comp_PLV([Y(1,:); Y(2,random_trials_indices)]);
    end
    
    %% Compute p-value
    % alpha=1e-6;
    % pct = quantile(PLV_permut,1-alpha);test statistic
    % plv_edge_decision=PLV_hat>pct;
    [f,x] = ecdf(PLV_permut);
    pvalue=1-max(f(x<=plv));
elseif d>2
    num_edges=d*(d-1)/2;
    
    plv =nan(num_edges,1);
    all_edges=nan(num_edges,2);
    pvalue=nan(num_edges,1);%rayleigh p values for plv
    inc=1;
    for j=1:d
        for k = (j+1):d
            plv(inc) = comp_PLV(Y([j,k],:));
            pvalue(inc) = circ_rtest(Y(j,:)-Y(k,:));
            all_edges(inc,:)=[j,k];
            inc=inc+1;
        end
    end
end



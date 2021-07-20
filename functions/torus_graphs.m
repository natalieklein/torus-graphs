function [TG,edges,phi_hat,inference]=torus_graphs(X,alpha_level,di,sel_mode,print_out,varargin)

%input description:
%-----------------
% X, d by N, matrix of angle vector of observations over repeated trials

% alpha_level, scalar value, probability of falsely placing an edge

% di, 1 by m, electrode grouping by m regions. d=sum(di)

% sel_mode is a boolean, 1 by 3, indicating which TG submodel to fit

% output description:
%-----------------
% TG,"graph with undirected edges" object, see https://www.mathworks.com/help/matlab/ref/graph.html}

% edges: structure with information about graph edges:
%   egdes.all_possible, num_edges by 2, array of all possible edges between node pairs
%   egdes.p_vals, num_edges by 1, edge p-values corresponding to egdes.all_possible
%   edges.active, active edges pval<alpha_level
%   edges.cond_coupling_coeff, transformation of the parameters that falls
%        between [0,1], phase difference model with uniform margins, see Section 4.4 in the paper.

% phi_hat, num_param by 1, vector of estimated natural parameters

% inference: structure with information to obtain p-values and perform
% inference on edges
%   Sigma_hat, num_param by num_param, covariance for phi
%   t_hat, num_edges by 1, test statistic for each edge group
%   dofs, num_edges by 1, number of degrees of freedom for each edge group
%%
[d,N]=size(X);

% Default input arguments
if ~exist('alpha_level','var') || isempty(alpha_level)
    alpha_level=0.05;
end

if ~exist('di','var') || isempty(di)
    di=ones(1,d); %one electrode within each brain region
end

if ~exist('print_out','var') || isempty(print_out)
    print_out=false;
end

if ~exist('sel_mode','var') || isempty(sel_mode)
    sel_mode=[true true true];
end

assert(sum(sel_mode)>0,'sel_mode, 000 not allowed')
[all_edges,params4edges,dofs_edges]=param_indexing(di,sel_mode);
%% Estimation
% [Gamma_hat_sub, H_hat_sub, Gamma_sub, H_sub,inds_from_full]=dx_suf_stat(X,sel_mode);
[Gamma_hat_sub, H_hat_sub,inds_from_full,V_hat]=dx_suf_stat2(X,sel_mode,print_out);

phi_hat_sub=Gamma_hat_sub\H_hat_sub;
phi_hat=zeros(2*d^2,1);
phi_hat(inds_from_full)=phi_hat_sub;
% [J,dJ]=SM_cost_fun(phi_hat_sub,[],Gamma_hat_sub, H_hat_sub);
% figure; stem(phi_hat)
%% Asymptotic covariance matrix and inference
% num_param_sub=size(phi_hat_sub,1);
% V_sum=zeros(num_param_sub,num_param_sub);
% for m=1:n
%     V_sum=V_sum+(Gamma_sub(:,:,m)*phi_hat_sub -H_sub(:,m))*(Gamma_sub(:,:,m)*phi_hat_sub -H_sub(:,m))';
% end
% V_hat=(1/n)*V_sum;
% Sigma_phi_hat=(1/n)*inv(Gamma_hat)*V_hat*inv(Gamma_hat);
Sigma_phi_hat_sub=(1/N)*(Gamma_hat_sub\V_hat)/Gamma_hat_sub;
Sigma_hat=zeros(2*d^2,2*d^2);
Sigma_hat(inds_from_full,inds_from_full)=Sigma_phi_hat_sub;
% figure;imagesc(Sigma_phi_hat);caxis([-.01 .01])

%%
if sel_mode(2) || sel_mode(3)
    num_edges=size(all_edges,1);
    t_hat=nan(num_edges,1);
    t_hat_pos=nan(num_edges,1);
    t_hat_neg=nan(num_edges,1);
    cond_coupling_squared=nan(num_edges,1);
    for inc=1:num_edges
        inds_edge=params4edges{inc};
        %     t_hat(inc)=phi_hat(inds_edge)'*inv(Sigma_phi_hat(inds_edge,inds_edge))*phi_hat(inds_edge);
        t_hat(inc)=(phi_hat(inds_edge)'/Sigma_hat(inds_edge,inds_edge))*phi_hat(inds_edge);
        if ~sel_mode(1) && sel_mode(2) && ~sel_mode(3) % phase difference model with uniform margins.
            cond_coupling_squared(inc)=phi_hat(inds_edge)'*phi_hat(inds_edge);
        end
        phi_hat_this_edge=phi_hat(inds_edge);
        Sigma_hat_this_edge=Sigma_hat(inds_edge,inds_edge);
        
        %------ Warning, separate pos/neg edges are only defined for
        %individual edges here, so it won't work for groupd of edges------
        % to do: fix it so that if di ~=[] then this part is not returned
        % or think of a way to define this
        if sel_mode(2) && sel_mode(3) % compute pos and neg
            t_hat_pos(inc)=(phi_hat_this_edge(1:2)'/Sigma_hat_this_edge(1:2,1:2))*phi_hat_this_edge(1:2);
            t_hat_neg(inc)=(phi_hat_this_edge(3:4)'/Sigma_hat_this_edge(3:4,3:4))*phi_hat_this_edge(3:4);
        elseif sel_mode(2) % only have pos
            t_hat_pos(inc)=(phi_hat_this_edge(1:2)'/Sigma_hat_this_edge(1:2,1:2))*phi_hat_this_edge(1:2);
            t_hat_neg(inc)=nan;
        else % sel_mode(3) % only have neg
            t_hat_pos(inc)=nan;
            t_hat_neg(inc)=(phi_hat_this_edge(1:2)'/Sigma_hat_this_edge(1:2,1:2))*phi_hat_this_edge(1:2);
        end
        
%         w_hat(inc)=phi_hat(inds_edge)'*phi_hat(inds_edge);
%         w_hat_ste(inc)=2*sqrt(phi_hat(inds_edge)'*Sigma_hat(inds_edge,inds_edge)*phi_hat(inds_edge));
    end
    all_pvals = 1 - chi2cdf(t_hat,dofs_edges);
    edges_on=all_edges(all_pvals<=alpha_level,:);
    
    %% plot network
    TG = graph(edges_on(:,1), edges_on(:,2));
% 
else
    t_hat=[];
    all_pvals=[];
    edges_on=[];
    TG=[];
end

edges.all_possible=all_edges;
edges.p_vals=all_pvals;
edges.active=edges_on;
if ~sel_mode(1) && sel_mode(2) && ~sel_mode(3) % phase difference model with uniform margins.
    edges.cond_coupling_coeff=besseli(1,sqrt(cond_coupling_squared))./besseli(0,sqrt(cond_coupling_squared));
end

inference.Sigma_hat=Sigma_hat;
inference.t_hat=t_hat;
inference.dofs_edges=dofs_edges;



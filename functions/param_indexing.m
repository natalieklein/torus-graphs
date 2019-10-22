function [all_edges,params4edges,dofs4edges]=param_indexing(di,sel_mode,varargin)

if ~exist('sel_mode','var')
    sel_mode=[true true true];
end

d=sum(di);
num_groups=length(di);
num_edges=num_groups*(num_groups-1)/2;

%% Edge indexing
phi_inds_interactions=(2*d+1):2*d^2;
phi_inds_interactions=reshape(phi_inds_interactions,[],4);
alpha_mat_inds=pairwise_list2dxd_mat(phi_inds_interactions(:,1));
beta_mat_inds=pairwise_list2dxd_mat(phi_inds_interactions(:,2));
gamma_mat_inds=pairwise_list2dxd_mat(phi_inds_interactions(:,3));
delta_mat_inds=pairwise_list2dxd_mat(phi_inds_interactions(:,4));
%%
all_edges=zeros(num_edges,2);
% inds_list_mat=cell(nump);
params4edges=cell(num_edges,1); % parameters corresponding to each edge
dofs4edges=nan(num_edges,1); % degrees of freedom corresponding to each edge
group_lines=[0 cumsum(di)];
row_inc=1;
for i=1:num_groups
    inds_rows=(group_lines(i)+1):group_lines(i+1);
    
    for j=(i+1):num_groups
        inds_cols=(group_lines(j)+1):group_lines(j+1);
        
        all_edges(row_inc,:)=[i j];
        
        if sel_mode(2) || sel_mode(3) 
            if sel_mode(2) && sel_mode(3) % _ 1 1 include both positive and negative correlations
                edge_inds=[alpha_mat_inds(inds_rows,inds_cols), beta_mat_inds(inds_rows,inds_cols), gamma_mat_inds(inds_rows,inds_cols), delta_mat_inds(inds_rows,inds_cols)];
            elseif sel_mode(2) &&  ~sel_mode(3) % _ 1 0 include only positive correlations
                edge_inds=[alpha_mat_inds(inds_rows,inds_cols), beta_mat_inds(inds_rows,inds_cols)];
            else % _ 0 1 include only negative correlations
                edge_inds=[gamma_mat_inds(inds_rows,inds_cols), delta_mat_inds(inds_rows,inds_cols)];
            end
            params4edges{row_inc}=edge_inds(:);
            dofs4edges(row_inc)=length(edge_inds(:));
        else
            params4edges={};
            dofs4edges=[];
        end
        row_inc=row_inc+1;
    end
end
% inds_list_mat=pairwise_list2dxd_mat(params4edges);
% params4edges2=dxd_mat2pairwise_list(inds_list_mat)

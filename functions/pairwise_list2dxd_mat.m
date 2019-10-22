function dxd_mat=pairwise_list2dxd_mat(pairwise_list)
% pairwise_list is nump x 1
% nump is d*(d-1)/2
nump=size(pairwise_list,1);
d=0.5*(1 + sqrt(1 + 8*nump));

dxdMat_indexing=reshape(1:(d*d),d,d);
inds_edges_pairs_in_dxdMat=triu(dxdMat_indexing,1)';
inds_edges_pairs_in_dxdMat=inds_edges_pairs_in_dxdMat(tril(true(d),-1));

if iscell(pairwise_list)
    A=cell(d);
    A=cellfun(@(a) 0, A, 'Uniform', false);
    A(inds_edges_pairs_in_dxdMat)=pairwise_list;
    A=cellfun(@(a, b) a + b, A, A', 'Uniform', false);
else
    A=zeros(d);
    A(inds_edges_pairs_in_dxdMat)=pairwise_list;
    A=A+A';
end

dxd_mat=A;
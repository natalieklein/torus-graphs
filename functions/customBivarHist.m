function [surfHandle,Nh,Ch]=customBivarHist(X,nbins)
% X is 2 by N
xEdge=linspace(-pi,pi,nbins);
yEdge=linspace(-pi,pi,nbins);
h=histogram2(X(1,:)',X(2,:)',xEdge,yEdge,'normalization','probability','FaceColor','flat','edgecolor','none','ShowEmptyBins','on');
bin_vals=h.Values(:);
% color axis displays 95% of the density range and thus excludes color outliers
caxis([prctile(bin_vals,2.5), prctile(bin_vals,97.5)])
[Nh,Ch]=hist3(X','edges',{xEdge,yEdge});
axis([-pi pi -pi pi])
view(0,90);axis square
surfHandle = get(gca, 'child');
colormap('parula');
function [Xdif, Xsum]=pairwise(X)
%input:
% X is chanels by observation trials (d x n)
%output
% Xdif and Xsum, paiwise differences and sums (nump x n)
% rows are ordered as: 1-2, 1-3, 1-4, ..., (d-1)-d.

[d,n]=size(X);

nump=d*(d-1)/2;
Xdif=nan(nump,n);
Xsum=nan(nump,n);

row_inc=1;
for i=1:d
    Xi=X(i,:);
    for j=(i+1):d
        Xj=X(j,:);
        Xdif(row_inc,:)=wrapToPi(Xi-Xj);
        Xsum(row_inc,:)=wrapToPi(Xi+Xj);
        row_inc=row_inc+1;
    end
end


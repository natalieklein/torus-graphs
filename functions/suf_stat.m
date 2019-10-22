function [Sc,Ss, Salpha, Sbeta, Sgamma, Sdelta]=suf_stat(X)
%input:
% X is chanels by observation trials (d x n)
%output
% Sufficient statistics:
% Sc and Ss are (d x n) corresponding to the individual concentrations,
% Salpha, Sbeta, Sgamma, and Sdelta are (nump x n), number of pairs by
% trials
% rows are ordered as: 1-2, 1-3, 1-4, ..., (d-1)-d.
[Xdif, Xsum]=pairwise(X);
Sc=cos(X);
Ss=sin(X);
Salpha=cos(Xdif);
Sbeta=sin(Xdif);
Sgamma=cos(Xsum);
Sdelta=sin(Xsum);
function [Gamma_hat, H_hat,inds_from_full,V_hat]=dx_suf_stat2(X,sel_mode,print_out,varargin)
% Computes the derivative with respect to x terms of the sufficient statistics
%input:
% X is chanels by observation trials (d x n)
%output:
% Gamma_hat and H_hat are the terms in the Score Matching objective
% function

if ~exist('sel_mode','var')
    sel_mode=[true true true];
end

[d,n]=size(X);
nump=d*(d-1)/2;

%% Getting number of parameters and indexing
H=[];
inds_from_full=[];
[Sc,Ss, Salpha, Sbeta, Sgamma, Sdelta]=suf_stat(X(:,1));
if sel_mode(1) % whether to fit marginal concentrations
    H=[H; Sc; Ss];
    inds_from_full=[inds_from_full 1:(2*d)];
end
if sel_mode(2) % whether to fit positive correlations
    H=[H; 2*Salpha; 2*Sbeta];
    inds_from_full=[inds_from_full (2*d+(1:2*nump)) ];
end
if sel_mode(3) % whether to fit negative correlations
    H=[H; 2*Sgamma; 2*Sdelta];
    inds_from_full=[inds_from_full (2*d+2*nump+(1:2*nump)) ];
end
num_param=size(H,1);
%%
Gamma_sum=zeros(num_param,num_param);
H_sum=zeros(num_param,1);

for m=1:n
    if print_out
        m
    end
    [Sc,Ss, Salpha, Sbeta, Sgamma, Sdelta]=suf_stat(X(:,m));

    Dx_Sc=-diag(Ss);
    Dx_Ss=diag(Sc);
    Dx_Salpha=zeros(nump,d);
    Dx_Sbeta=zeros(nump,d);
    Dx_Sgamma=zeros(nump,d);
    Dx_Sdelta=zeros(nump,d);
    
    row_ind=1;
    for i=1:d
        for j=(i+1):d
            col_inds=[i,j];
            Dx_Salpha(row_ind,col_inds)=[-1 1]*Sbeta(row_ind);
            Dx_Sbeta(row_ind,col_inds)=[1 -1]*Salpha(row_ind);
            Dx_Sgamma(row_ind,col_inds)=[-1 -1]*Sdelta(row_ind);
            Dx_Sdelta(row_ind,col_inds)=[1 1]*Sgamma(row_ind);
            row_ind=row_ind+1;
        end
    end
    
    Dx_temp=[];
    H=[];
    if sel_mode(1) % whether to fit marginal concentrations
        H=[H; Sc; Ss];
        Dx_temp=[Dx_temp; Dx_Sc;Dx_Ss];
    end
    if sel_mode(2) % whether to fit positive correlations
        H=[H; 2*Salpha; 2*Sbeta];
        Dx_temp=[Dx_temp; Dx_Salpha;Dx_Sbeta];
    end
    if sel_mode(3) % whether to fit negative correlations
        H=[H; 2*Sgamma; 2*Sdelta];
        Dx_temp=[Dx_temp; Dx_Sgamma;Dx_Sdelta];
    end
    H_sum=H_sum + H;
    Gamma_sum=Gamma_sum + Dx_temp*Dx_temp';
end

H_hat=(1/n)*H_sum;
Gamma_hat=(1/n)*Gamma_sum;
phi_hat=Gamma_hat\H_hat;

V_sum=zeros(num_param,num_param);
for m=1:n
    if print_out
        m
    end
    [Sc,Ss, Salpha, Sbeta, Sgamma, Sdelta]=suf_stat(X(:,m));

    Dx_Sc=-diag(Ss);
    Dx_Ss=diag(Sc);
    Dx_Salpha=zeros(nump,d);
    Dx_Sbeta=zeros(nump,d);
    Dx_Sgamma=zeros(nump,d);
    Dx_Sdelta=zeros(nump,d);
    
    row_ind=1;
    for i=1:d
        for j=(i+1):d
            col_inds=[i,j];
            Dx_Salpha(row_ind,col_inds)=[-1 1]*Sbeta(row_ind);
            Dx_Sbeta(row_ind,col_inds)=[1 -1]*Salpha(row_ind);
            Dx_Sgamma(row_ind,col_inds)=[-1 -1]*Sdelta(row_ind);
            Dx_Sdelta(row_ind,col_inds)=[1 1]*Sgamma(row_ind);
            row_ind=row_ind+1;
        end
    end
    
    Dx_temp=[];
    H=[];
    if sel_mode(1) % whether to fit marginal concentrations
        H=[H; Sc; Ss];
        Dx_temp=[Dx_temp; Dx_Sc;Dx_Ss];
    end
    if sel_mode(2) % whether to fit positive correlations
        H=[H; 2*Salpha; 2*Sbeta];
        Dx_temp=[Dx_temp; Dx_Salpha;Dx_Sbeta];
    end
    if sel_mode(3) % whether to fit negative correlations
        H=[H; 2*Sgamma; 2*Sdelta];
        Dx_temp=[Dx_temp; Dx_Sgamma;Dx_Sdelta];
    end
    Gamma=Dx_temp*Dx_temp';
    
    V_sum=V_sum+(Gamma*phi_hat -H)*(Gamma*phi_hat -H)';
end
V_hat=(1/n)*V_sum;
function [Xsim,Xoriginaldata]=simAltTG(simname,N)
% simname is 'A' or 'B'
% N is desired sample size
load('data/anglebeta.mat')

ca3inds = find(ismember(subregions,'CA3'));
%
if strcmp(simname,'A')
    %% Simulation A: Spatial autoregressive 5-dim, compare to CA3 channels
    dinds = ca3inds;
    Xoriginaldata = wrapTo2Pi(anglebeta(dinds,:))';
    [n,d] = size(Xoriginaldata);
    n_noise = 15;
    n_samp = N; 

    %offset = pi/50;
    offset = [pi/100 pi/100 pi/100 pi/100 pi/100];
    noisekappa = 40;
    
    Xsim = zeros(n_samp,d+1);
    Xsim(:,1) = wrapTo2Pi(circ_vmrnd(0,0.01,n_samp));
    for i=2:d+1
        exnoise = wrapTo2Pi(circ_vmrnd(0,0.1,n_samp));
        exnoise(randsample(n_samp,n_samp-n_noise)) = 0;
        Xsim(:,i) = wrapTo2Pi(Xsim(:,i-1) + offset(i-1) + circ_vmrnd(0,noisekappa,n_samp) + exnoise);
        %X(:,i) = wrapTo2Pi(X(:,i-1) + offset + pi/25*(rand(1,1)-0.5) + circ_vmrnd(0,noisekappa,n_samp) + exnoise);
    end
    
    Xsim = Xsim(:,2:end);
elseif strcmp(simname,'Aalt')
    dinds = ca3inds;
    Xoriginaldata = wrapTo2Pi(anglebeta(dinds,:))';
    [n,d] = size(Xoriginaldata);
    n_noise = 25;
    n_samp = N; 
    noisekappa = 50;
    
    offsets = pi/50*[-1,1;-1,1;-1,1;-1,1;-1,1;-1,1;-1,1;-1,1;-1,1];
    %offsets = pi/5*[1,1;1,1;1,1;1,1;1,1];
    
    % independent theta
    %theta = zeros(n_samp,4);
    %for i=1:4
    %    theta(:,i) = wrapTo2Pi(circ_vmrnd(0,0.01,n_samp));
    %end
    
    % correlated theta
    theta = zeros(n_samp,d+1);
    theta(:,1) = wrapTo2Pi(circ_vmrnd(0,0.01,n_samp));
    for i=2:d+1
        theta(:,i) = theta(:,i-1) + pi/100 + circ_vmrnd(0,noisekappa,n_samp);
    end
    
    exnoise = wrapTo2Pi(circ_vmrnd(0,0.01,n_samp));
    exnoise(randsample(n_samp,n_samp-n_noise)) = 0;
    
    Xsim = zeros(n_samp,d+2);
    Xsim(:,1) = theta(:,1) + offsets(1,1) + circ_vmrnd(0,noisekappa,n_samp);
    for i=2:d+1
        Xtmp = theta(:,i-1) + theta(:,i) + offsets(i-1,2) + offsets(i,1)  + circ_vmrnd(0,noisekappa,n_samp);
        %Xtmp = theta(:,i-1) + offsets(i-1,2) + circ_vmrnd(0,noisekappa,n_samp) + exnoise;
        %Xtmp = Xtmp + theta(:,i) + offsets(i,1) + circ_vmrnd(0,noisekappa,n_samp) + exnoise;
        Xsim(:,i) = Xtmp;
    end
    Xsim(:,d+2) = theta(:,d+1) + offsets(d+1,2) + circ_vmrnd(0,noisekappa,n_samp);
    
    Xsim = wrapTo2Pi(bsxfun(@plus,Xsim(:,2:d+1),exnoise));
    
elseif strcmp(simname,'B')
    %% Simulation B: confounder
    % Selected these based on a 24dim SM graph for having partial correlations
    dinds = [13,16,23];
    Xoriginaldata = wrapTo2Pi(anglebeta(dinds,:))'; % 13=DG,16=Sub,23=PFC
    [n,d] = size(Xoriginaldata);
    
    n_noise = 150;
    n_samp = N; 
    Xsim = zeros(N,d);
    
    offset1 = pi/6;
    offset2 = pi/100;
    noisekappa = 3;

    % Generate Sub
    Xsim(:,2) = wrapTo2Pi(circ_vmrnd(0,0.01,n_samp));
    % Generate DG from Sub
    exnoise = circ_vmrnd(0,0.1,n_samp);
    exnoise(randsample(n_samp,n_samp-n_noise)) = 0;
    Xsim(:,1) = wrapTo2Pi(Xsim(:,2) + offset1 + circ_vmrnd(0,noisekappa,n_samp) + exnoise);
    % Generate PFC from Sub
    exnoise = circ_vmrnd(0,0.1,n_samp);
    exnoise(randsample(n_samp,n_samp-n_noise)) = 0;
    Xsim(:,3) = wrapTo2Pi(Xsim(:,2) + offset2 + circ_vmrnd(0,noisekappa,n_samp) + exnoise);
end

Xsim=wrapToPi(Xsim);
Xoriginaldata=wrapToPi(Xoriginaldata);


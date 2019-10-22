function S = sampleGibbs_tutorial(d,phi,options)
% Samples using Gibbs sampling
% Params
%   options : struct with (nsamp,burnin,thin) options
% Returns
%   S       : d by nsamp matrix of samples

% d = this.dim;
totalSamples = options.burnin + options.nsamp*options.thin;
S = zeros(d,totalSamples);

% Initialize uniformly
x = 2*pi*rand(1,d);


%             function [mu,kappa,Cm,Sm,Cp,Sp] = phiToParams(this,phi)
% d = this.dim;
Cm = zeros(d,d);
Sm = zeros(d,d);
Cp = zeros(d,d);
Sp = zeros(d,d);

Ne=d*(d-1)/2;% number of edges
inds=tril(ones(d,d),-1)==1;
phiCminds=(2*d+1):(2*d+Ne);
phiSminds=(2*d+Ne+1):(2*d+2*Ne);
phiCpinds=(2*d+2*Ne+1):(2*d+3*Ne);
phiSpinds=(2*d+3*Ne+1):(2*d^2);
Cm(inds) = phi(phiCminds);
Cm = Cm + Cm';
Sm(inds) = phi(phiSminds);
Sm = Sm + Sm';
Cp(inds) = phi(phiCpinds);
Cp = Cp + Cp';
Sp(inds) = phi(phiSpinds);
Sp = Sp + Sp';
cosmu = phi(1:d);
sinmu = phi((d+1):2*d);
kappa = sqrt(cosmu.^2+sinmu.^2);
mu = wrapTo2Pi(atan2(sinmu,cosmu));
%             end

% Extract parameters
tmu = mu;
tkappa = kappa;
tCp = Cp;
tSp = Sp;
tCm = Cm;
tSm = Sm;

for samp=1:totalSamples
    for k=1:d % get sample of k conditional on rest
        % old version; changing signs to match document version of theorem
        %deltavals = [-tmu(k),-x(1:d~=k),-x(1:d~=k)-pi/2,x(1:d~=k),x(1:d~=k)-pi/2];
        %Avals = [tkappa(k),tCm(k,1:d~=k),-tSm(k,1:d~=k),tCp(k,1:d~=k),tSp(k,1:d~=k)];
        %[A,delta] = harmonicAddition(Avals,deltavals);
        %x(k) = wrapToPi(circ_vmrnd(-delta,A,1));
        
        % version as of 2/28: gives same as above
        %                     deltavals = [tmu(k),x(1:d~=k),x(1:d~=k)-pi/2,-x(1:d~=k),-x(1:d~=k)+pi/2];
        %                     Avals = [tkappa(k),tCm(k,1:d~=k),tSm(k,1:d~=k),tCp(k,1:d~=k),tSp(k,1:d~=k)];
        %                     [A,delta] = harmonicAddition(Avals,deltavals);
        %                     %x(k) = wrapToPi(circ_vmrnd(delta,A,1));
        %                     x(k) = circ_vmrnd(delta,A,1);
        
        % fixing on 3/16
        smdelta = zeros(1,d-1);
        smdelta(1:d<k) = x(1:d<k)-pi/2;
        smdelta(find(1:d>k)-1) = x(1:d>k)+pi/2;
        deltavals = [tmu(k),x(1:d~=k),smdelta,-x(1:d~=k),-x(1:d~=k)+pi/2];
        Avals = [tkappa(k),tCm(k,1:d~=k),tSm(k,1:d~=k),tCp(k,1:d~=k),tSp(k,1:d~=k)];
        [A,delta] = harmonicAddition(Avals,deltavals);
        %x(k) = wrapToPi(circ_vmrnd(delta,A,1));
        x(k) = circ_vmrnd(delta,A,1);
    end
    S(:,samp) = x;
end
S = S(:,options.burnin+1:options.thin:end);
end
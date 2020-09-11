function [nNf] = TLModel_GOK(time,temp,kparams);

Ma = 3600*24*365*1e6;
% Extract parameters
ddot = kparams(1)/(1000*365*24*3600);                  % Gy/ka to Gy/seconds
D0 = kparams(2);
a  = kparams(3);
Et = kparams(4);
% Eu = kparams(5);
s = 10.^kparams(6);
b = kparams(7);
rhop = 10.^kparams(8);


% Define constants
kb = 8.617343e-5; alpha = 1; Hs = 3e15; %s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
nrp = 100;
% nEb = 50;
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma; %time step in Myr is converted into second.

% Define rprime range and tau athermic
rprime = linspace(0.01,2.5,nrp); %create vector of rprime distances
pr = 3.*rprime'.^2.*exp(-rprime'.^3.); %calc p(r') eq 3 in Kars et al 2008
npr = sum(pr);

inv_tauath = (Hs*exp(-(rhop.^-(1./3)).*rprime))'; %combine eq 1 and 3 from Kars et al 2008, convert to Ma

% Create variables for lili2013 BTM
% Eb = zeros(1,nEb)';
% Eb = linspace(0.01,Et+0.01,nEb)';
% pEb = exp(-Eb'/Eu);
% npEb = sum(pEb);
T = temp+273.15;

% computes nN for a Tt path
nN = zeros(nrp,nstep);
nNf = zeros(1,nstep);
for i = 2:nstep
	inv_tauth = (s*exp(-(Et)./(kb.*T(i-1)))*ones(1,nrp))';
    xkd = -a*magic_ratio*(1-nN(:,i-1)).^(a-1)-b*(nN(:,i-1)).^(b-1).*inv_tauth-inv_tauath;
    xk = magic_ratio*(1-nN(:,i-1)).^a-(nN(:,i-1)).^b.*inv_tauth-inv_tauath.*nN(:,i-1);
    nN(:,i) = nN(:,i-1)+xk.*dt./(1-dt.*xkd);
% 	alpha = magic_ratio+inv_tauth+inv_tauath;
%     nN(:,:,i) = (nN(:,:,i-1)+magic_ratio.*dt).*(1./(1+dt.*alpha));
	nNf(i) = nN(:,i)'*pr;
end
nNf = nNf./npr;

% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)

clear all; close all; clc;
tic
%% Load TL parameters
filename = 'Summary_Parameter_GOK_MBTP9';
rawdata=xlsread(filename);

%  pool=parpool(10);

kparams=rawdata(:,1:2:18); skparams=rawdata(:,2:2:18);


load(['MBTP9_Path1_nN_obs']);
nN_mean = nN_obs; sigmanN=0.2*nN_mean;

%% Parameters for Random Path1s Path1
timemax=1.0;        % Time in Ma
timemin=0;          % Time in Ma
nstep=2000;         % this is to reinterpolate the TtPath1 for the nN model
niter=1000;         % number of random Path1s, fixed to equal the number of random realisations

Dtime   =timemax/(nstep-1);
timeM   =0:Dtime:timemax;


parfor i=1:niter
residuals=zeros(size(kparams,1),1);
nNf=zeros(nstep,size(kparams,1),niter);
Predicted_nN=zeros(niter,size(kparams,1));
nN=zeros(size(kparams,1),1);
residuals=zeros(size(kparams,1),1);
nNf=zeros(nstep,size(kparams,1),niter);
Predicted_nN=zeros(niter,size(kparams,1));
nN=zeros(size(kparams,1),1);


% Generating Random oscillation Path1 
meanT   = 50*rand;
amp     = 50*rand;                            % in degree C
    period_Mil  = 1/1000*[25.77,41,100];
    period  = period_Mil(randperm(3,1));
% period=25.77/1000;
tempM1=meanT+amp*cos(2*pi/period*timeM);
tempM =fliplr(tempM1);


   
plot(timeM, tempM); hold on
    xlim([max(timeM)-100/1000,max(timeM)]);                  % Plotting of last 200 ka
% ylim([-10,30]);
    xlabel('time (ka)');
    ylabel('temp (^oC)');
    set(gca,'FontSize',15);

%% Calculate nN for each Tt Path1s
 
    for k=1:size(kparams,1)
        nNf(:,k,i) = TLModel_GOK_MCerror(timeM,tempM,kparams(k,:),skparams(k,:));
        nN(k)=nN_mean(k)+2*sigmanN(k)*(rand-0.5);
        residuals(k)=((nN(k)/sigmanN(k))*0.5.*log(nN(k)/nNf(end,k,i))).^2;
        Predicted_nN(:,k)=nNf(end,k,i);    
    end
    sum_residual=sum(residuals)/length(nN);
    
    misfit(i,1)=sum_residual;
    TPath1(i,:)=tempM;
    nNPred(i,:) = Predicted_nN(i,:);
    
    Tmean(i,1)=meanT;
    Tamp(i,1)=amp;
    Tperiod(i,1)=period;
    
    

fprintf('Path1%i   \n ',i);
end
%%
Out_misfit = misfit; 
Out_TPath1 = TPath1;
Out_nNPred = nNPred;
Out_Tmean = Tmean;
Out_Tamp = Tamp;
Out_Tperiod = Tperiod;


% Sorting the data
disp('sorting table')
[sortedmisOUT,IX] = sort(Out_misfit(:,1));
sortedTt = Out_TPath1(IX,:);
sortedPrednN = Out_nNPred(IX,:);
sortedTmean=Out_Tmean(IX,1);
sortedTamp=Out_Tamp(IX,1);
sortedTperiod=Out_Tperiod(IX,1);

% Saving the data
save('MBTP9_Path1_Tt.mat', 'sortedTt', '-v7.3');
save 'MBTP9_Path1_misfit.mat' sortedmisOUT
save 'MBTP9_Path1_PrednN.mat' sortedPrednN
save 'MBTP9_Path1_time.mat' timeM
save 'MBTP9_Path1_Tmean.mat' sortedTmean
save 'MBTP9_Path1_Tamp.mat' sortedTamp
save 'MBTP9_Path1_Tperiod.mat' sortedTperiod
% delete(pool);
toc

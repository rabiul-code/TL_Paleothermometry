% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)

clear all; close all; clc;

%% Load TL parameters
load('MBTP9_TL_params.mat');
TLtemp = 205:10:295;            % TL signals that need to invert
for i=1:length(TLtemp)
[sig(i),idx] = find(TL_params(:,1)==TLtemp(i));
end
Im      =TL_params(sig,2);
kparams =TL_params(sig,3:end);

%% Parameters for Random oscillation path
timeM    = linspace(0, 1, 2000);      % In ka
meanT   = [15 30];
amp     = [10 20];
period  = 1/1000*[10 20];                         % in ka

% temp = meanT(4)*ones(length(timeM),1);
% for k=1:size(kparams,1)                        % number of cooling Paths
%     nNf = TLModel_GOK(timeM,temp,kparams(k,:));
%     Out_nN_isoT(k,:)=nNf;
%     Out_nNmax_isoT(k,:)=nNf(end);
% end


tempM(:,1)=meanT(1)+amp(1)*sin(2*pi/period(1)*timeM);
tempM(:,2)=meanT(2)+amp(1)*sin(2*pi/period(1)*timeM);
tempM(:,3)=meanT(1)+amp(2)*sin(2*pi/period(1)*timeM);
tempM(:,4)=meanT(1)+amp(1)*sin(2*pi/period(2)*timeM);
for i=1:4
    
    for k=1:size(kparams,1)                        % number of cooling Paths
        nNf1 = TLModel_GOK(timeM,tempM(:,i),kparams(k,:));
        Out_nN(i,k,:)=nNf1;
        Out_nNmax(k,i)=nNf1(end);
    end
end

f2 = figure(2); axis square; box on; hold on
% plot(TLtemp,Out_nNmax_isoT,'o','MarkerSize', 12);
plot(TLtemp,Out_nNmax,'o-','MarkerSize', 8,'Linewidth', 1);
xlim([200 300]);
ax=gca;
set(gca,'YLim',[1e-3 5],'YTick',[10e-5 10e-4 10e-3 10e-2 1]);
ax.YAxis.TickLabelFormat = '10^-%d';

xlabel('TL temperature (^oC)');
ylabel('n');
set(gca,'FontSize',20);
set(gca, 'Yscale', 'log')
% LEG=legend({'T_{mean} = 0 ^oC', 'T_{mean} = 15 ^oC','T_{mean} = 30 ^oC'}, 'Location', 'SouthEast');
LEG=legend({'T_{amp} = 5 ^oC','T_{amp} = 10 ^oC','T_{amp} = 20 ^oC'}, 'Location', 'SouthEast');
% LEG=legend({'P = 1 ka', 'P = 10 ka','P = 100 ka'}, 'Location', 'SouthEast');

legend boxoff;
LEG.FontSize = 20;
% A=title(['T_{amp} = 10 ^oC;     P = 10 ka']);
A=title(['T_{mean} = 15 ^oC;     P = 10 ka']);
% A=title(['T_{amp} = 10 ^oC;     T_{mean} = 15 ^oC']);

A.FontSize = 15;

% print(f2, 'MBTP9_Sensitivity_Amp', '-dpdf', '-r300');

% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)

clear all; close all; clc;
set(gca, 'FontName', 'Times')

%% Load TL parameters
load('MBTP9_TL_params.mat');
TLtemp = 205:20:295;            % TL signals that need to invert
for i=1:length(TLtemp)
[sig(i),idx] = find(TL_params(:,1)==TLtemp(i));
end
Im      =TL_params(sig,2);
kparams =TL_params(sig,3:end);

%% Temp define: define what kind of temp field
timeM    = linspace(0, 1, 10000)';      

% Oscillating field
n           =3;
amp     = 10;
meanT       =30;
period  = 100/1000;      % in Ma
tempM_isoT  =meanT.*ones(length(timeM),1);
tempM_period   =meanT+amp*sin(2*pi/period*timeM);

%% Plotting of temp field
f1=figure(1); axis square; box on; hold on 
pbaspect([9 9 1]);
col = [0.3 0.3 0.3; 0.4660    0.6740    0.1880; 0.8500    0.3250    0.0980];
for i=1:length(meanT)
plot(timeM,tempM_period(:,i), 'col', col(i,:), 'LineWidth', 1.0);
plot(timeM,tempM_isoT(:,i), '--','col', col(i,:),'LineWidth', 1.5);
end
xlim([max(timeM)-300/1000,max(timeM)]);
% ylim([-10,30]);
xlabel('time (ka)');
ylabel('temp (^oC)');
set(gca,'FontSize',27);
set(gca,'XTickLabel', {'300','200','100','0'});
print('Input thermal field_100ka', '-dpdf', '-r300');
%% Calculating n/N
for k=1:size(kparams,1)                        % number of cooling Paths
    nN_isoT(:,k) = TLModel_GOK(timeM,tempM_isoT,kparams(k,:));
    nN_period(:,k) = TLModel_GOK(timeM,tempM_period,kparams(k,:));
end

f2=figure(2); axis square; box on; hold on
pbaspect([1 2 1]);
% cc =colormap(othercolor('Reds7',size(kparams,1))); 
cc = [0.3 0.3 0.3; 0.4660    0.6740    0.1880; 0.8500    0.3250    0.0980];

for k=1:size(kparams,1)
plot(timeM,nN_period(:,k),'col',cc(n,:),'LineWidth', 1.0);
plot(timeM,nN_isoT(:,k),'--','col',cc(n,:),'LineWidth', 1.0);
end
set(gca,'Yscale','log');
xlim([max(timeM)-300/1000,max(timeM)]);
ylim([8*1e-4 2]);
xlabel('time (ka)');
ylabel('n');
set(gca,'FontSize',17);
set(gca,'XTickLabel', {'300','200','100','0'});



% print(f2,'MBTP9_nN_evolution_100ka_3', '-dpdf', '-r300');

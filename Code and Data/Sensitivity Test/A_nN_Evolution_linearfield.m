% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)

clear all; close all; clc;
set(gca, 'FontName', 'Times')

%% Load TL parameters
load('MBTP9_TL_params.mat');
TLtemp = 205:10:295;            % TL signals that need to invert
for i=1:length(TLtemp)
[sig(i),idx] = find(TL_params(:,1)==TLtemp(i));
end
Im      =TL_params(sig,2);
kparams =TL_params(sig,3:end);

%% Temp define: define what kind of temp field
timeM    = linspace(0, 1, 2000)';      

% Isothermal field
meanT       =20;
n           =3;
tempM_isoT  =meanT.*ones(length(timeM),1);

% Cooling field
cut_time = 1-100/1000;
temp_max = 20;
temp_min = 0;
x1=[min(timeM),cut_time,max(timeM)];
y1=[temp_max,temp_max,temp_min];
tempM_cooling=interp1(x1,y1,timeM);
% 
% Warming field
cut_time = 1-100/1000;
temp_max = 20;
temp_min = 0;
x2=[min(timeM),cut_time,max(timeM)];
y2=[temp_min,temp_min,temp_max];
tempM_warming=interp1(x2,y2,timeM);


%% Plotting of temp field
f1=figure(1); axis square; box on; hold on 
pbaspect([9 9 1]);
col = [0.3 0.3 0.3; 0.4660    0.6740    0.1880; 0.8500    0.3250    0.0980];
for i=1:length(meanT)
plot(timeM,tempM_isoT(:,i),'color',[0.5 0.5 0.5],'LineWidth', 3);
end
xlim([max(timeM)-200/1000,max(timeM)]);
ylim([-10,30]);
xlabel('time (ka)');
ylabel('temp (^oC)');
set(gca,'FontSize',27);
set(gca,'XTickLabel', {'200','100','0'});
print('Input thermal field_100ka', '-dpdf', '-r300');
%% Calculating n/N
for k=1:size(kparams,1)                        % number of isoT Paths
    nN_isoT(:,k) = TLModel_GOK(timeM,tempM_isoT,kparams(k,:));
end

f2=figure(2); axis square; box on; hold on
pbaspect([1 2 1]);
% cc =colormap(othercolor('Reds7',size(kparams,1))); 
% cc = [0.3 0.3 0.3; 0.4660    0.6740    0.1880; 0.8500    0.3250    0.0980];
cc=othercolor('Reds9',12);
for k=1:size(kparams,1)
% plot(timeM,nN_period(:,k),'col',cc(n,:),'LineWidth', 1.0);
plot(timeM,nN_isoT(:,k),'color',cc(k+1,:),'LineWidth', 2.5);
end
set(gca,'Yscale','log');
xlim([max(timeM)-200/1000,max(timeM)]);
ylim([1e-3 2]);
xlabel('time (ka)');
ylabel('n');
set(gca,'FontSize',17);
set(gca,'XTickLabel', {'200','100','0'});



% print(f1,'MBTP9_Path_isoT', '-dpdf', '-r300');
% print(f2,'MBTP_nN_isoT', '-dpdf', '-r300');


% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)

clear all; close all; clc;

%% Load TL parameters
filename = 'Summary_Parameter_GOK_MBTP9';
rawdata=xlsread(filename);

% pool=parpool(10);

kparams=rawdata(:,1:2:18); skparams=rawdata(:,2:2:18);

%% Parameters for Random Path1
%% Defining Parameters
timeM    = linspace(0, 1, 4000);      % In ka
meanT   = 10;
amp     = 10;                            % in degree C
period  = 25.77/1000;                             % in Ma
tempM1=meanT+amp*cos(2*pi/period*timeM);
tempM =fliplr(tempM1);
Path1(:,1)=timeM;
Path1(:,2)=tempM;
save(['MBTP9_Path1.mat'],'Path1');

%% Main Program
for k=1:size(kparams,1)                        % number of cooling Path1s
    nNf = TLModel_GOK(timeM,tempM,kparams(k,:));
    Out_nN_time(k,:)=nNf;
    nN_obs(k,:)=nNf(end);
end
save (['MBTP9_Path1_nN_obs.mat'], 'nN_obs');
%% Plotting
f1=figure(1); axis square; box on; hold on 
% pbaspect([16 9 1]);
plot1=plot(timeM,tempM,'color',[0.5 0.5 0.5], 'LineWidth', 3);
plot1.Color(4) = 0.8;
xlim([max(timeM)-100/1000,max(timeM)]);                  % Plotting of last 200 ka
ylim([-10,40]);
set(gca,'XTickLabel', {'100','80','60','40','20','0'});
xlabel('time (ka)');
ylabel('Temperature (^oC)');
set(gca,'FontSize',20);

f2=figure(2); axis square; box on; hold on
cc=othercolor('Greys9',8);
for i=1:size(kparams,1)
    plot(timeM,Out_nN_time(i,:),'LineWidth', 1.0,'color',cc(i+2,:),'LineWidth',3.0)
    set(gca,'Yscale','log');
    xlim([max(timeM)-100/1000,max(timeM)]);
    ylim([1e-3 2]);
    set(gca,'XTickLabel', {'100','80','60','40','20','0'});
    xlabel('time (ka)');
    ylabel('n');
    set(gca,'FontSize',20);
end
% 
% print(f1,'MBTP9_Forward_Path1', '-dpdf', '-r300');
% print(f2,'MBTP9_Forward_nN_Path1', '-dpdf', '-r300');
% 


    
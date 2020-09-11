% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)

clear all; close all; clc;

%% Load TL parameters
filename = 'Summary_Parameter_GOK_MBTP9';
rawdata=xlsread(filename);

% pool=parpool(10);

kparams=rawdata(:,1:2:18); skparams=rawdata(:,2:2:18);

%% Parameters for Random Path1
%% Defining Parameters
timemax=1.0;        % Time in Ma
timemin=0;          % Time in Ma
nstep=2000;         % this is to reinterpolate the TtPath1 for the nN model

Dtime   =timemax/(nstep-1);
timeM   =0:Dtime:timemax;

% temp define from delta 18O data
data = xlsread('Greenland_delta18O');
time =data(:,1)/10^6;            % time in Ma
detal18O=data(:,2);

% time_scaled    = linspace(0, 1, 2000);      % In ka
detal18O_f =fliplr(interp1(time,detal18O,timeM));

%scaling
[val,idx]=min(abs(timeM-0.98)); %base point at 20 ka
max_temp=detal18O_f(end);
min_temp=detal18O_f(idx);

detal18O_norm=(-min_temp+detal18O_f)./(max_temp-min_temp);

meanT=0;
amp=10;


tempM=meanT+detal18O_norm.*amp;

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
plot1=plot(timeM,tempM,'color',[0.3 0.3 0.3], 'LineWidth', 3);
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

% print(f1,'MBTP9_Forward_Path1', '-dpdf', '-r300');
% print(f2,'MBTP9_Forward_nN_Path1', '-dpdf', '-r300');



    
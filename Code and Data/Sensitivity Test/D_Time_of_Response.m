% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)

close all; clear all;
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');

%% Load TL parameters
load('MBTP9_TL_params.mat');
TLtemp = 205:10:295;            % TL signals that need to invert
for i=1:length(TLtemp)
[sig(i),idx] = find(TL_params(:,1)==TLtemp(i));
end
Im      =TL_params(sig,2);
kparams =TL_params(sig,3:end);

%% Define Tt Paths
timeM   =linspace(0,1,1000)';
max_region = 150/1000;                              % 100 ka
step       = 1/1000;                               % 10 ka
t_step  =max(timeM)-max_region:step:max(timeM);
temp_b4   =0;
temp_step =30;

for j=1:length(t_step)
    for i=1:length(timeM)
        if timeM(i)<t_step(j)
            tempM(i,j)=temp_b4;
        else
            tempM(i,j)=temp_step;
        end
    end
end


%% Calculate nN
nN_present =zeros(size(kparams,1),length(t_step));
for j=1:length(t_step)                             % number of paths
    for k=1:size(kparams,1)                        % number of signals
    nN_time = TLModel_GOK(timeM,tempM(:,j),kparams(k,:));
    nN_present(k,j)=nN_time(end);
    end
end
f1=figure(1); axis square; box on; hold on
cc=othercolor('Reds9',12);
for i=1:size(kparams,1)
% plot(t_step,nN_present(i,:),'o','MarkerSize',5,'MarkerFaceColor',cc(i+4,:),'MarkerEdgeColor','none');hold on
plot(t_step,nN_present(i,:),'color',cc(i+1,:),'LineWidth',3.0);
[val(i),idx(i)]=min(abs(1.1*nN_present(i,1)-nN_present(i,:)));

% [val(i),idx(i)]=min(abs(0.1*(nN_present(i,end)-nN_present(i,1))+nN_present(i,1)-nN_present(i,:)));
end

for i=1:size(kparams,1)+1
    plot(t_step(idx(i)),nN_present(i,idx(i)),'k*','MarkerSize',15);
    if idx(i+1)>idx(i)   
        break
    end
end    

xlim([0.85 1.0]);
ylim([0.001 9]);
set(gca,'XTickLabel', {'100','80','60','40','20','0'});
set(gca,'Yscale','log')
axes.XDir = 'reverse';
xlabel('t_{change} (ka)');
ylabel('Present day n');
set(gca,'FontSize',18);








% print(f1,'MBTP9_time_response_Oto30', '-dpdf', '-r300');
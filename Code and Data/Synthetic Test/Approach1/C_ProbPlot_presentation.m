% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)


tic

clear all;clc;clf,close all,
set(gcf,'units','centimeters')
set(gca,'FontSize',10)
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');

% Loading the file
filename='MBTP9_Path1';
load([filename '_Tt.mat']); Tt=sortedTt;
load([filename '_misfit.mat']); misOUT=sortedmisOUT;     
% load([filename '_nN.mat']); nNs=sortednNsave; 
load([filename '_PrednN.mat']); PrednNs=sortedPrednN; 
load([filename '_time.mat']); time=timeM;

load([filename '_Tmean.mat']); Tmean=sortedTmean;
load([filename '_Tamp.mat']); Tamp=sortedTamp;
load([filename '_Tperiod.mat']); Tperiod=sortedTperiod*1000;

load([filename '_nN_obs']);
nN = nN_obs; sigmanN=0.2*nN;

MTemp = [215 225 235 245];            % TL signals that need to invert

% Grid define
nAvM=1000;
[m,nt]=size(PrednNs);

% to compute the PDF for the likelhood, the Zt paths are reinterpolated onto a grid (Av_matrix)
time_max=1; time_min=0; dt=(time_max-time_min)/(nAvM-1);
T_max=500; T_min=-500; dT=(T_max-T_min)/(nAvM-1);
vec_time=time_min:dt:time_max;
Tvec=T_min:dT:T_max;

% Data selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = rand(m,1); 
prob=exp(-misOUT); scale=max(prob); % scaling for resampling the PDF
% scale=1;
test=prob/scale>R;
%test=scale/prob>R;

idefix=find(test);
idefix=idefix(end:-1:1); % goes from largest misfit to smallest
movea=length(idefix);

max_likelihood=max(prob(idefix)); min_likelihood=min(prob(idefix)); % bounds within the max and min you want the color scheme to be
index_color=max(1,floor(63*(prob-min_likelihood)/(max_likelihood-min_likelihood))+1); % scaling the value between 1 and 64

MTEMP=ones(movea,1)*MTemp(1:nt);
INDEX=index_color(idefix)*ones(1,nt);
PrednNs=PrednNs(idefix,:);

Av_matrix=zeros(nAvM);
% add accepted model to a matrix to compute the PDF, this is meand on a rejection algorithm
for k = idefix';
	vec_T=interp1(time,Tt(k,:),vec_time,'linear');
    Tpath=(0:nAvM-1)*nAvM+round((vec_T-T_min)/dT)+1;
    Tpath(Tpath<=0)=[];
	Av_matrix(Tpath)=Av_matrix(Tpath)+1;
end
X=cumsum(Av_matrix/movea); % for computing CIs and median


% Plot figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map=colormap(parula); % creates a table with 64 rows and three columns for RGB
col=round(100*(max(misOUT(idefix)):(min(misOUT)-max(misOUT(idefix)))/6:min(misOUT)))/100; % create colorbar for misfit

% Figure 1: Predicted vs. observed nN
f1=figure(1); axis square; box on; hold on
xlabel('TL temperature (^oC)');
ylabel('n');
c1=colorbar('yticklabel',col);
set(get(c1,'title'),'string','misfit');
colormap(parula);

median_path=scatter(MTEMP(:),PrednNs(:),60,INDEX(:),'filled');   % predicted nNs
actual=errorbar(MTemp(1:nt),nN(1:nt),sigmanN(1:nt),'ko','MarkerSize', 20);   % observed nNs
legend([actual,median_path],'Observed','Predicted', 'Location', 'NorthWest'); legend boxoff;

xlim([200 250]);
ylim([1e-3 10]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
set(gca,'Yscale','log');


% Figure 2: Spaghetti plot of selected paths
f2=figure(2); axis square; box on; hold on
xlabel('Time (Ma)');
ylabel('Temp (^oC)');
c1=colorbar('yticklabel',col);
colormap(parula);
set(get(c1,'title'),'string','misfit');
P = plot(time,Tt(idefix,:)); % 'spaghetti' plot
for curve = 1:movea
	set(P(curve),'Color',map(index_color(idefix(curve)),:))
end

xlim([time_max-0.1 time_max]);
ylim([-20 40])
set(gca,'FontSize',20);
ax.LineWidth = 2.0;

% Figure 3: Probability Density Plot
f3=figure(3); axis square; box on; hold on
xlabel('Time (ka)');
ylabel('Temperature (^oC)');
axis([time_max-0.05 time_max -20 50],'square')
c2=colorbar;colormap(parula);
set(get(c2,'title'),'string','PDF');
hold on
contourf(vec_time,Tvec,Av_matrix./movea,100,'edgecolor','none'); % PDF
shading flat
twosig=contour(vec_time,Tvec,X,[0.05,0.95],'k','LineWidth',0.5); % 95 CI
onesig=contour(vec_time,Tvec,X,[0.32,0.68],'w','LineWidth',1.5); % 68 CI
median_path=contour(vec_time,Tvec,X,1,'r','LineWidth',2); % median
load(['MBTP9_Path1.mat']);
actual=plot(Path1(:,1), Path1(:,2), 'r:','LineWidth',3.0)

caxis([0 0.1]); 
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
set(gca,'XTickLabel', {'50','40','30','20','10','0'});
% set(gca,'xaxisLocation','top');
% set(gca,'yaxisLocation','right');
% colorbar('westoutside')
% set(gca,'caxisLocation','west');


Tmean_final=Tmean(idefix,1);            n_Tmean=round(max(Tmean_final)-min(Tmean_final));
Tamp_final=Tamp(idefix,1);              n_Tamp=round(max(Tamp_final)-min(Tamp_final));


% Mean temp histogram
f4=figure(4); axis square; box on; hold on
P1=histogram(Tmean_final,n_Tmean,'FaceColor',[0.2 0.2 0.2]);
[hi cx]=hist(Tmean_final,n_Tmean);
AA=median(Tmean_final);
P2=bar(AA,max(hi));
P3=bar(10,max(hi));
xlabel('Mean Temp (^oC)');
ylabel('Frequency');
xlim([-20 30]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
legend([P1, P2, P3], 'Histogram','Predicted median', 'Actual')


% Amplitude histogram
f5=figure(5); axis square; box on; hold on
P1=histogram(Tamp_final,n_Tamp,'FaceColor',[0.2 0.2 0.2]);
[hi cx]=hist(Tamp_final,n_Tamp);
AA=median(Tamp_final);
P2=bar(AA,max(hi));
P3=bar(10,max(hi));
xlabel('Amplitude (^oC)');
ylabel('Frequency');
% xlim([0 40]);
legend([P1, P2, P3], 'Histogram','Predicted median', 'Actual')
set(gca,'FontSize',20);
ax.LineWidth = 2.0;

% Period histogram
Tperiod_final=Tperiod(idefix,1);        n_Tperiod=round(max(Tperiod_final)-min(Tperiod_final));
f6=figure(6); axis square; box on; hold on
histogram(Tperiod_final,n_Tperiod,'FaceColor',[0.2 0.2 0.2]);
[hi cx]=hist(Tperiod_final,n_Tperiod);
AA=median(Tperiod_final);
% bar(AA,max(hi));
xlabel('Period (ka)');
ylabel('Frequency');
% xlim([0 40]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;


% print(f1,'MBTP9_Path1_nN_prediction', '-dpdf', '-r300');
% print(f3,'MBTP9_Path1_probPlot', '-dpng', '-r300');
% print(f4,'MBTP9_Path1_Tmean', '-dpdf', '-r300');
% print(f5,'MBTP9_Path1_Tamp', '-dpdf', '-r300');
% print(f6,'MBTP9_Path1_Tperiod', '-dpdf', '-r300');

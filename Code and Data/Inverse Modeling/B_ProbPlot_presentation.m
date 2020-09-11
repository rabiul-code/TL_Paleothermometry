% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)

tic

clear all;clc;clf,close all,
set(gcf,'units','centimeters')
set(gca,'FontSize',10)
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');

% Loading the file
filename='MBTP1';
load([filename '_Tt.mat']); Tt=sortedTt;
load([filename '_misfit.mat']); misOUT=sortedmisOUT;
% load([filename '_nN.mat']); nNs=sortednNsave; 
load([filename '_PrednN.mat']); PrednNs=sortedPrednN; 
load([filename '_time.mat']); time=timeM;

load([filename '_Tbase.mat']); Tbase=sortedTbase;
load([filename '_Tamp.mat']); Tamp=sortedTamp;

rawdata=xlsread('Summary_Parameter_GOK_MBTP1');
nN=rawdata(:,23); sigmanN=rawdata(:,24);
KarsnN=rawdata(:,25); sigmaKarsnN=rawdata(:,26);

MTemp=[215 225 235 245];

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
% add accepted model to a matrix to compute the PDF, this is based on a rejection algorithm
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
% colormap(parula);

median_path=scatter(MTEMP(:),PrednNs(:),60,INDEX(:),'filled');   % predicted nNs
actual=errorbar(MTemp(1:nt),nN(1:nt),sigmanN(1:nt),'ko','MarkerSize', 20);   % observed nNs
legend([actual,median_path],'Observed','Predicted', 'Location', 'NorthWest'); legend boxoff;

xlim([200 250]);
ylim([1e-2 10]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
set(gca,'Yscale','log');


% Figure 2: Spaghetti plot of selected paths
f2=figure(2); axis square; box on; hold on
xlabel('Time (Ma)');
ylabel('Temp (^oC)');
c1=colorbar('yticklabel',col);colormap(jet);
set(get(c1,'title'),'string','misfit');
P = plot(time,Tt(idefix,:)); % 'spaghetti' plot
for curve = 1:movea
	set(P(curve),'Color',map(index_color(idefix(curve)),:))
end

xlim([time_max-0.1 time_max]);
ylim([-20 20])
set(gca,'FontSize',20);
ax.LineWidth = 2.0;

% Figure 3: Probability Density Plot
f3=figure(3); axis square; box on; hold on
xlabel('Time (ka)');
ylabel('Temperature (^oC)');
axis([time_max-0.05 time_max -20 20],'square');
cc=othercolor('Reds9',100)
c2=colorbar;
set(get(c2,'title'),'string','PDF');
hold on
contourf(vec_time,Tvec,Av_matrix./movea,100,'edgecolor','none'); % PDF
shading flat
twosig=contour(vec_time,Tvec,X,[0.05,0.95],'k','LineWidth',0.5); % 95 CI
onesig=contour(vec_time,Tvec,X,[0.32,0.68],'w','LineWidth',1.5); % 68 CI
median_path=contour(vec_time,Tvec,X,1,'r','LineWidth',2); % median

caxis([0 0.1]); 
set(gca,'FontSize',20);
ax.LineWidth = 2.0;
set(gca,'XTickLabel', {'50','40','30','20','10','0'});
% set(gca,'xaxisLocation','top');
% set(gca,'yaxisLocation','right');
% colorbar('westoutside')
% set(gca,'caxisLocation','west');


Tbase_final=Tbase(idefix,1);            n_Tbase=round(max(Tbase_final)-min(Tbase_final));
Tamp_final=Tamp(idefix,1);              n_Tamp=round(max(Tamp_final)-min(Tamp_final));
Tsurf_final=Tbase_final+Tamp_final;     n_Tsurf=round(max(Tsurf_final)-min(Tsurf_final));


f4=figure(4); axis square; box on; hold on
histogram(Tbase_final,n_Tbase,'FaceColor',[0.4 0.4 0.4]);
[hi cx]=hist(Tbase_final,n_Tbase);
AA=median(Tbase_final);
bar(AA,max(hi));
xlabel('Temp at 20 ka (^oC)');
ylabel('Frequency');
xlim([-20 30]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;

f5=figure(5); axis square; box on; hold on
histogram(Tsurf_final,n_Tsurf,'FaceColor',[0.2 0.2 0.2]);
[hi cx]=hist(Tsurf_final,n_Tsurf);
BB=median(Tsurf_final);
bar(BB,max(hi));
xlabel('Temp at present (^oC)');
ylabel('Frequency');
xlim([-20 30]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;


f6=figure(6); axis square; box on; hold on
histogram(Tamp_final,n_Tamp,'FaceColor',[0.2 0.2 0.2]);
[hi cx]=hist(Tamp_final,n_Tamp);
AA=median(Tamp_final);
bar(AA,max(hi));
xlabel('Tamp (^oC)');
ylabel('Frequency');
% xlim([0 40]);
set(gca,'FontSize',20);
ax.LineWidth = 2.0;


f7=figure(7); axis square; box on; hold on
xlabel('Tbase (^oC)');
ylabel('Tamp (^oC)');
c1=colorbar('yticklabel',col);colormap(jet);
set(get(c1,'title'),'string','misfit');
for i=1:movea
    P(i) = plot(Tbase_final(i),Tamp_final(i),'o'); % 'spaghetti' plot
end

for curve = 1:movea
	set(P(curve),'Color',map(index_color(idefix(curve)),:))
end

xlim([-20 30]);
ylim([0 40])
set(gca,'FontSize',20);
ax.LineWidth = 2.0;

%%
% Median and one sigma path
outliers = find (median_path(2,:) > 100);
median_path(:,outliers)=[];
[~,idx]=find(median_path(1,:)==1);
median_path=median_path(:,1:idx);

outliers = find (onesig(2,:) > 100);
onesig(:,outliers)=[];

[~,idx]=find(onesig(1,:)==1);
onesig_minus=onesig(:,1:idx(1));
onesig_plus=onesig(:,1+idx(1):end);

% Present day temp
T_present_median=median_path(2,end);
T_present_plus=onesig_plus(2,end);
T_present_minus=onesig_minus(2,end);


% LGM temp
LGM1=1-20/1000;
LGM2=1-17/1000;
LGM3=1-15/1000;
LGM4=1-12/1000;

[~,idx]=min(abs(median_path(1,:)-LGM1)); T_LGM1_median=median_path(2,idx);
[~,idx]=min(abs(onesig_plus(1,:)-LGM1)); T_LGM1_plus=onesig_plus(2,idx);
[~,idx]=min(abs(onesig_minus(1,:)-LGM1)); T_LGM1_minus=onesig_minus(2,idx);

[~,idx]=min(abs(median_path(1,:)-LGM2)); T_LGM2_median=median_path(2,idx);
[~,idx]=min(abs(onesig_plus(1,:)-LGM2)); T_LGM2_plus=onesig_plus(2,idx);
[~,idx]=min(abs(onesig_minus(1,:)-LGM2)); T_LGM2_minus=onesig_minus(2,idx);

[~,idx]=min(abs(median_path(1,:)-LGM3)); T_LGM3_median=median_path(2,idx);
[~,idx]=min(abs(onesig_plus(1,:)-LGM3)); T_LGM3_plus=onesig_plus(2,idx);
[~,idx]=min(abs(onesig_minus(1,:)-LGM3)); T_LGM3_minus=onesig_minus(2,idx);

[~,idx]=min(abs(median_path(1,:)-LGM4)); T_LGM4_median=median_path(2,idx);
[~,idx]=min(abs(onesig_plus(1,:)-LGM4)); T_LGM4_plus=onesig_plus(2,idx);
[~,idx]=min(abs(onesig_minus(1,:)-LGM4)); T_LGM4_minus=onesig_minus(2,idx);


T_extract(1,1)=T_present_median;
T_extract(1,2)=T_present_plus;
T_extract(1,3)=T_present_minus;
T_extract(1,4)=T_LGM1_median;
T_extract(1,5)=T_LGM1_plus;
T_extract(1,6)=T_LGM1_minus;

T_extract(2,1)=T_present_median;
T_extract(2,2)=T_present_plus;
T_extract(2,3)=T_present_minus;
T_extract(2,4)=T_LGM2_median;
T_extract(2,5)=T_LGM2_plus;
T_extract(2,6)=T_LGM2_minus;


T_extract(3,1)=T_present_median;
T_extract(3,2)=T_present_plus;
T_extract(3,3)=T_present_minus;
T_extract(3,4)=T_LGM3_median;
T_extract(3,5)=T_LGM3_plus;
T_extract(3,6)=T_LGM3_minus;

T_extract(4,1)=T_present_median;
T_extract(4,2)=T_present_plus;
T_extract(4,3)=T_present_minus;
T_extract(4,4)=T_LGM4_median;
T_extract(4,5)=T_LGM4_plus;
T_extract(4,6)=T_LGM4_minus;


% print(f1,'MBTP1_nN_prediction', '-dpdf', '-r300');
% print(f3,'MBTP1_probPlot', '-dpdf', '-r300');
% print(f4,'MBTP1_Tbase', '-dpdf', '-r300');
% print(f5,'MBTP1_Tamp', '-dpdf', '-r300');
% % print(f6,'MBTP1_Tamp_vs_Tbase', '-dpdf', '-r300');

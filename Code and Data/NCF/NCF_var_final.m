% written by Rabiul Haque Biswas (biswasrabiul@gmail.com)

clear all; close all; clc;

file = xlsread('MBTP_NCF');
rawdata = file(2:end-1, 2:end);

for i=1:size(rawdata,2)       
    TL1=rawdata(:,1:2:end);
    TL2=rawdata(:,2:2:end);
    ratio =TL1./TL2;
end

temp=linspace(0,200,250)';
% plot(temp,TL2);

% Measured NCF
Temp_meas= 90:10:150;
ncf_meas=interp1(temp,ratio,Temp_meas);
for k=1:length(Temp_meas)
    for i=1:size(ncf_meas,2)/3
        ncf_meas_mean(k,i)=mean(ncf_meas(k,3*i-2:3*i));
        ncf_meas_std(k,i)=std(ncf_meas(k,3*i-2:3*i));
    end  
end

% Extrapolated NCF
Temp_extra=[205 215 225 235 245 255 265 275 285 295]';
Temp_plot=linspace(90,300,10);
for i=1:size(ncf_meas,2)
    P(i,:)=polyfit(Temp_meas',ncf_meas(:,i),1);
    ncf_extra(:,i)=P(i,1)*Temp_extra+P(i,2);
    ncf_plot(:,i)=P(i,1)*Temp_plot+P(i,2);
end
for k=1:length(Temp_extra)
    for i=1:size(ncf_meas,2)/3
        ncf_extra_mean(k,i)=mean(ncf_extra(k,3*i-2:3*i));
        ncf_extra_std(k,i)=std(ncf_extra(k,3*i-2:3*i));
    end  
end 



% Plotting TL1 and TL2 of 1st sample first disk
f1=figure(1);axis square; box on; hold on
plot(temp,TL1(:,1),'Linewidth',2.0);
plot(temp,TL2(:,1),'Linewidth',2.0);
xlim([0 400]);
% ylim([0 3]);
legend('Before Natural (TL1)','After Natural (TL2)');
legend boxoff;
xlabel('TL temperature (^oC)');
ylabel('TL Intensity (counts/s)');
set(gca,'FontSize',20);

% Plotting var-NCF of 1st sample 2nd disk
f2=figure(2); axis square; box on; hold on
col=[0    0.4470    0.7410; 0.8500    0.3250    0.0980];

% MBTP1 disc 2, 1st sample
P1=plot(Temp_meas,ncf_meas(:,3*1-1),'o','MarkerSize', 10, 'MarkerFaceColor',col(1,:),'MarkerEdgeColor','none');
plot(Temp_extra(2:5),ncf_extra(2:5,3*1-1),'o','MarkerSize', 10, 'MarkerFaceColor','none','MarkerEdgeColor',col(1,:));
plot(Temp_plot,ncf_plot(:,3*1-1),'Color',col(1,:));

% MBTP9 disc 2, 5th sample
P2=plot(Temp_meas,ncf_meas(:,3*5-1),'o','MarkerSize', 10, 'MarkerFaceColor',col(2,:),'MarkerEdgeColor','none');
plot(Temp_extra(2:5),ncf_extra(2:5,3*5-1),'o','MarkerSize', 10, 'MarkerFaceColor','none','MarkerEdgeColor',col(2,:));
plot(Temp_plot,ncf_plot(:,3*5-1),'Color',col(2,:));

% for i=1:2
%     P(i)=plot(Temp_meas,ncf_meas(:,3*i-1),'o','MarkerSize', 10, 'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none'); 
%     plot(Temp_extra(2:5),ncf_extra(2:5,3*i-1),'o','MarkerSize', 10, 'MarkerFaceColor','none','MarkerEdgeColor',col(i,:)); 
%     plot(Temp_plot,ncf_plot(:,3*i-1),'Color',col(i,:));
% end
xlim([90 260]);
ylim([0 3]);
legend([P1 P2],'MBTP1','MBTP9');
% legend boxoff;
xlabel('TL temperature (^oC)');
ylabel('var-NCF');
set(gca,'FontSize',20);


print(f1,'MBTP1_NCF', '-dpdf', '-r300');
print(f2,'MBTP1 and 9 var_NCF', '-dpdf', '-r300');

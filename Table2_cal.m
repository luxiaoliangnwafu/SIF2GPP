% Table 2

PAR1 = 600;
Tair1 =20;
%% GPPsif GPPec  half hourly  NDVI<0.6 

figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
%
% ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
% set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

DateStrings = {'2021-02-14'};
date1 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');
DateStrings = {'2021-05-30'};
date2 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');


axes(ha(1));

x_input_temp = mean([SIF_GPP_dataset1.GPP(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset2.GPP(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset3.GPP(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset4.GPP(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset5.GPP(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset6.GPP(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )],2,'omitnan');
                                    
y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 )],2,'omitnan');                                   
                                        
% y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime < t) SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset2.Datetime < t) SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset3.Datetime < t)...
%     SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset4.Datetime < t) SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset5.Datetime < t) SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset6.Datetime < t) ],2,'omitnan');

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input = SIF_GPP_dataset1.PAR_dataset0(SIF_GPP_dataset1.Datetime < date1 &  SIF_GPP_dataset1.Datetime < date2 );
x_input = t1;
y_input = t2;

y=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
x=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(x,y,360,color_bar,'.');
%set(get(get(g1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(h1, 'CLim', [min(color_input) max(color_input)]);
h1= colorbar('eastoutside');
colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
% h1.Limits = [0 	1000];
% h1.Ticks = [0:200:1000];
% plot(x_input,y_input,'ko',...
%     'MarkerSize',6,...
%     'MarkerEdgeColor',1/255.*[150 150 150],...
%     'MarkerFaceColor',1/255.*[150 150 150]);

% if max(x_input) > max(y_input)
%     xlim([0 1.1*max(x_input)])
%     ylim([0 1.1*max(x_input)])
% else
%     xlim([0 1.1*max(y_input) ])
%     ylim([0 1.1*max(y_input) ])
% end

% xlim([0 40])
% set(gca,'XTick',[0: 10 :40],'ycolor','k')
% ylim([0 40])
% set(gca,'YTick',[0: 10 :40],'ycolor','k')

xlabel(['GPP_S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
hold on
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
plot(xrange, xrange,'r-.','LineWidth',2) % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
plot(xrange,ypred,'k-','LineWidth',2)
boundedline(xrange, ypred, delta,'alpha');
% plot(xrange,delta_t,'r--','LineWidth',1)
% plot(xrange,ypred,'k-','LineWidth',1.5)
hold off
lgd = legend('Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',22);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(x_input)-min(x_input));
title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.48*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
%set(gca,'LooseInset',get(gca,'TightInset'));
box on

temp1 = {num2str(R22,'%6.2f');num2str(RMSE2,'%6.2f');num2str(rrmse,'%6.2f');num2str(size(x_input,1),'%6.0f')};

%% GPPsif GPPec  half hourly  NDVI>0.6 

figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
%
% ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
% set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);


DateStrings = {'2021-02-14'};
date1 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');
DateStrings = {'2021-05-30'};
date2 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');


axes(ha(1));

x_input_temp = mean([SIF_GPP_dataset1.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset2.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset3.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset4.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset5.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset6.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )],2,'omitnan');
                                    
y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )...
                                        SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 )],2,'omitnan');                                   
                                        
% y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime < t) SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset2.Datetime < t) SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset3.Datetime < t)...
%     SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset4.Datetime < t) SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset5.Datetime < t) SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset6.Datetime < t) ],2,'omitnan');

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input = SIF_GPP_dataset1.PAR_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 );
x_input = t1;
y_input = t2;

y=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
x=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(x,y,360,color_bar,'.');
%set(get(get(g1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(h1, 'CLim', [min(color_input) max(color_input)]);
h1= colorbar('eastoutside');
colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';


xlabel(['GPP_S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
hold on
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
plot(xrange, xrange,'r-.','LineWidth',2) % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
plot(xrange,ypred,'k-','LineWidth',2)
boundedline(xrange, ypred, delta,'alpha');
% plot(xrange,delta_t,'r--','LineWidth',1)
% plot(xrange,ypred,'k-','LineWidth',1.5)
hold off
lgd = legend('Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',22);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(x_input)-min(x_input));
x_input(1:8) = [];
title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.48*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
%set(gca,'LooseInset',get(gca,'TightInset'));
box on

temp2 = {num2str(R22,'%6.2f');num2str(RMSE2,'%6.2f');num2str(rrmse,'%6.2f');num2str(size(x_input,1),'%6.0f')};

%% GPPsif GPPec  half hourly  NDVI>0.6 PAR<600

figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
%
% ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
% set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

DateStrings = {'2021-02-14'};
date1 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');
DateStrings = {'2021-05-30'};
date2 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');


axes(ha(1));

x_input_temp = mean([SIF_GPP_dataset1.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset2.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset3.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset4.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset5.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset6.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)],2,'omitnan');
                                    
y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)...
                                        SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1)],2,'omitnan');                                   
                                        
% y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime < t) SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset2.Datetime < t) SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset3.Datetime < t)...
%     SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset4.Datetime < t) SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset5.Datetime < t) SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset6.Datetime < t) ],2,'omitnan');

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input = SIF_GPP_dataset1.PAR_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 <= PAR1);
x_input = t1;
y_input = t2;

y=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
x=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(x,y,360,color_bar,'.');
%set(get(get(g1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(h1, 'CLim', [min(color_input) max(color_input)]);
h1= colorbar('eastoutside');
colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
h1.Limits = [0 	1000];
h1.Ticks = [0:200:1000];

xlabel(['GPP_S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
hold on
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
plot(xrange, xrange,'r-.','LineWidth',2) % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
plot(xrange,ypred,'k-','LineWidth',2)
boundedline(xrange, ypred, delta,'alpha');
% plot(xrange,delta_t,'r--','LineWidth',1)
% plot(xrange,ypred,'k-','LineWidth',1.5)
hold off
lgd = legend('Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',22);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(x_input)-min(x_input));
title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.48*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
%set(gca,'LooseInset',get(gca,'TightInset'));
box on

temp3 = {num2str(R22,'%6.2f');num2str(RMSE2,'%6.2f');num2str(rrmse,'%6.2f');num2str(size(x_input,1),'%6.0f')};



%% GPPsif GPPec  half hourly  NDVI>0.6 PAR>600

figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
%
% ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
% set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);


DateStrings = {'2021-02-14'};
date1 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');
DateStrings = {'2021-05-30'};
date2 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');


axes(ha(1));

x_input_temp = mean([SIF_GPP_dataset1.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset2.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset3.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset4.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset5.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset6.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)],2,'omitnan');
                                    
y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)...
                                        SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1)],2,'omitnan');                                   
   
t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input = SIF_GPP_dataset1.PAR_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.PAR_dataset0 > PAR1);
x_input = t1;
y_input = t2;

y=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
x=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(x,y,360,color_bar,'.');
%set(get(get(g1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(h1, 'CLim', [min(color_input) max(color_input)]);
h1= colorbar('eastoutside');
colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
% h1.Limits = [0 	1000];
% h1.Ticks = [0:200:1000];
% plot(x_input,y_input,'ko',...
%     'MarkerSize',6,...
%     'MarkerEdgeColor',1/255.*[150 150 150],...
%     'MarkerFaceColor',1/255.*[150 150 150]);

% if max(x_input) > max(y_input)
%     xlim([0 1.1*max(x_input)])
%     ylim([0 1.1*max(x_input)])
% else
%     xlim([0 1.1*max(y_input) ])
%     ylim([0 1.1*max(y_input) ])
% end

% xlim([0 40])
% set(gca,'XTick',[0: 10 :40],'ycolor','k')
% ylim([0 40])
% set(gca,'YTick',[0: 10 :40],'ycolor','k')

xlabel(['GPP_S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
hold on
model_temp = fitlm([x_input],y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
plot(xrange, xrange,'r-.','LineWidth',2) % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
plot(xrange,ypred,'k-','LineWidth',2)
boundedline(xrange, ypred, delta,'alpha');
% plot(xrange,delta_t,'r--','LineWidth',1)
% plot(xrange,ypred,'k-','LineWidth',1.5)
hold off
lgd = legend('Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',22);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(x_input)-min(x_input));
title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.48*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
%set(gca,'LooseInset',get(gca,'TightInset'));
box on
% txt2 = { 'Half hourly'} ;
% text(0.4*max(xlim),1.08*max(ylim),txt2,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
% txt = {'d'};
% text(0.95*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',32)
%%%print('D:\Dropbox\Fluorescence_Liu_10\figure\Figure5\Figure5_YL_GPP_halfhourly_winter','-dtiff','-r300')
temp4 = {num2str(R22,'%6.2f');num2str(RMSE2,'%6.2f');num2str(rrmse,'%6.2f');num2str(size(x_input,1),'%6.0f')};


%% GPPsif GPPec  half hourly  NDVI>0.6 Tair<20

figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
%
% ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
% set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);


DateStrings = {'2021-02-14'};
date1 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');
DateStrings = {'2021-05-30'};
date2 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');



axes(ha(1));

x_input_temp = mean([SIF_GPP_dataset1.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset2.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset3.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset4.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset5.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset6.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)],2,'omitnan');
                                    
y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)...
                                        SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1)],2,'omitnan');                                   
                                        
% y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime < t) SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset2.Datetime < t) SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset3.Datetime < t)...
%     SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset4.Datetime < t) SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset5.Datetime < t) SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset6.Datetime < t) ],2,'omitnan');

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input = SIF_GPP_dataset1.PAR_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 <= Tair1);
x_input = t1;
y_input = t2;

y=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
x=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(x,y,360,color_bar,'.');
%set(get(get(g1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(h1, 'CLim', [min(color_input) max(color_input)]);
h1= colorbar('eastoutside');
colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
% h1.Limits = [0 	1000];
% h1.Ticks = [0:200:1000];
% plot(x_input,y_input,'ko',...
%     'MarkerSize',6,...
%     'MarkerEdgeColor',1/255.*[150 150 150],...
%     'MarkerFaceColor',1/255.*[150 150 150]);

% if max(x_input) > max(y_input)
%     xlim([0 1.1*max(x_input)])
%     ylim([0 1.1*max(x_input)])
% else
%     xlim([0 1.1*max(y_input) ])
%     ylim([0 1.1*max(y_input) ])
% end

% xlim([0 40])
% set(gca,'XTick',[0: 10 :40],'ycolor','k')
% ylim([0 40])
% set(gca,'YTick',[0: 10 :40],'ycolor','k')

xlabel(['GPP_S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
hold on
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
plot(xrange, xrange,'r-.','LineWidth',2) % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
plot(xrange,ypred,'k-','LineWidth',2)
boundedline(xrange, ypred, delta,'alpha');
% plot(xrange,delta_t,'r--','LineWidth',1)
% plot(xrange,ypred,'k-','LineWidth',1.5)
hold off
lgd = legend('Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',22);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(x_input)-min(x_input));
title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.48*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
%set(gca,'LooseInset',get(gca,'TightInset'));
box on
% txt2 = { 'Half hourly'} ;
% text(0.4*max(xlim),1.08*max(ylim),txt2,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
% txt = {'d'};
% text(0.95*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',32)
%%%print('D:\Dropbox\Fluorescence_Liu_10\figure\Figure5\Figure5_YL_GPP_halfhourly_winter','-dtiff','-r300')

temp5 = {num2str(R22,'%6.2f');num2str(RMSE2,'%6.2f');num2str(rrmse,'%6.2f');num2str(size(x_input,1),'%6.0f')};




%% GPPsif GPPec  half hourly  NDVI>0.6 Tair>20
%axes(ha(4));
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
%
% ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
% set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

DateStrings = {'2021-02-14'};
date1 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');
DateStrings = {'2021-05-30'};
date2 = datetime(DateStrings,'InputFormat','yyyy-MM-dd');

axes(ha(1));

x_input_temp = mean([SIF_GPP_dataset1.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset2.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset3.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset4.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset5.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset6.GPP(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)],2,'omitnan');
                                    
y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)...
                                        SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1)],2,'omitnan');                                   
                                        
% y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0(SIF_GPP_dataset1.Datetime < t) SIF_GPP_dataset2.GPP_dataset0(SIF_GPP_dataset2.Datetime < t) SIF_GPP_dataset3.GPP_dataset0(SIF_GPP_dataset3.Datetime < t)...
%     SIF_GPP_dataset4.GPP_dataset0(SIF_GPP_dataset4.Datetime < t) SIF_GPP_dataset5.GPP_dataset0(SIF_GPP_dataset5.Datetime < t) SIF_GPP_dataset6.GPP_dataset0(SIF_GPP_dataset6.Datetime < t) ],2,'omitnan');

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input = SIF_GPP_dataset1.PAR_dataset0(SIF_GPP_dataset1.Datetime > date1 &  SIF_GPP_dataset1.Datetime < date2 & SIF_GPP_dataset1.Temperature_dataset0 > Tair1);
x_input = t1;
y_input = t2;

y=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
x=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
g1 = scatter(x,y,360,color_bar,'.');
%set(get(get(g1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(h1, 'CLim', [min(color_input) max(color_input)]);
h1= colorbar('eastoutside');
colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
% h1.Limits = [0 	1000];
% h1.Ticks = [0:200:1000];
% plot(x_input,y_input,'ko',...
%     'MarkerSize',6,...
%     'MarkerEdgeColor',1/255.*[150 150 150],...
%     'MarkerFaceColor',1/255.*[150 150 150]);

% if max(x_input) > max(y_input)
%     xlim([0 1.1*max(x_input)])
%     ylim([0 1.1*max(x_input)])
% else
%     xlim([0 1.1*max(y_input) ])
%     ylim([0 1.1*max(y_input) ])
% end

% xlim([0 40])
% set(gca,'XTick',[0: 10 :40],'ycolor','k')
% ylim([0 40])
% set(gca,'YTick',[0: 10 :40],'ycolor','k')

xlabel(['GPP_S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
hold on
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
plot(xrange, xrange,'r-.','LineWidth',2) % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
plot(xrange,ypred,'k-','LineWidth',2)
boundedline(xrange, ypred, delta,'alpha');
% plot(xrange,delta_t,'r--','LineWidth',1)
% plot(xrange,ypred,'k-','LineWidth',1.5)
hold off
lgd = legend('Half hourly','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',22);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(x_input)-min(x_input));
title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ¡À ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.48*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
%set(gca,'LooseInset',get(gca,'TightInset'));
box on

temp6 = {num2str(R22,'%6.2f');num2str(RMSE2,'%6.2f');num2str(rrmse,'%6.2f');num2str(size(x_input,1),'%6.0f')};

table2= [temp1 temp2 temp3 temp4 temp5 temp6];

close all









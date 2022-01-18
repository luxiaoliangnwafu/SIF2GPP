clc;clear;
load RSE_submit.mat
%% GPP calculation
GPP_datetime = datetime(GPP_dataset.Year,1,GPP_dataset.DoY)+hours(GPP_dataset.Hour);
for ii=1:6    
    % run six observed points in Yangling sites    
    eval(['SIF_dataset0=','SIF_dataset',num2str(ii),';'])
    
    % clear value before the growing wheat
    index  = find (SIF_dataset0.SVD_O2A>0);
    SIF_dataset0(1: index(1)-1,:)=[];
    clear index
    
    % find the paired GPP data
    [~, idx] = intersect(GPP_datetime, SIF_dataset0.Datetime);
    GPP_dataset0 = GPP_dataset.GPP_DT(idx);
    Temperature_dataset0 = GPP_dataset.Tair(idx);
    PAR_dataset0 = GPP_dataset.PPFD(idx);
    rH =str2double(GPP_dataset.Rh(idx))./100;
    VPD = (1-rH).*(0.6107*exp(17.27.*Temperature_dataset0./(Temperature_dataset0+237.3)));
    Ca = GPP_dataset.CO2_mean(idx);
    vwc = GPP_dataset.VWC_mean(idx)./100;
    
    %  clear value equal or below zero
    index  = find (SIF_dataset0.SVD_O2A<=0 | PAR_dataset0<=20 | GPP_dataset0<=0);
    SIF_dataset0(index,:)=[];
    GPP_dataset0(index)=[];
    PAR_dataset0(index)=[];
    Temperature_dataset0(index)=[];
    Ca(index)=[];
    VPD(index)=[];
    vwc(index)=[];
    clear index
    
    % outliers detect and replace
    t = SIF_dataset0.Datetime;
    A = GPP_dataset0;
    [B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(6),'SamplePoints',t);
    GPP_dataset0 = B;
    clear A B t
    
    t = SIF_dataset0.Datetime;
    A = PAR_dataset0;
    [B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(6),'SamplePoints',t);
    PAR_dataset0 = B;
    clear A B t
    
    t = SIF_dataset0.Datetime;
    A = Temperature_dataset0;
    [B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(12),'SamplePoints',t);
    Temperature_dataset0 = B;
    clear A B t
    
    t = SIF_dataset0.Datetime;
    A = VPD;
    [B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(6),'SamplePoints',t);
    VPD = B;
    clear A B t
    
    t = SIF_dataset0.Datetime;
    A = Ca;
    [B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(6),'SamplePoints',t);
    Ca = B;
    clear A B t
    
    t = SIF_dataset0.Datetime;
    A = vwc;
    [B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(6),'SamplePoints',t);
    vwc = B;
    clear A B t
    
    % phi PSII Bacur et al., 2019 JGR
    % ETR_j    
    a0=0.9;%BACOUR 2019 Eq5
    b0=0.611;
    c0=12.132;
    d0=-0.086;
    Cab=50;
    alpha = a0.*(1-b0.*exp(-Cab./c0+d0));  %BACOUR 2019 Eq5
    
    t = SIF_dataset0.Datetime;
    A = SIF_dataset0.fPAR;
    [B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(12),'SamplePoints',t);
    %fAPAR_dataset0 = B;
    SIF_dataset0.fPAR =B;
    clear A B t
    fesc = 1./(SIF_dataset0.fPAR./SIF_dataset0.NIR_v_O2A./SIF_dataset0.NDVI);
    
    t = SIF_dataset0.Datetime;
    A = fesc;
    [B,~,~,~,~] = filloutliers(A,'Linear','movmedian',hours(6),'SamplePoints',t);
    fesc = B;
    clear A B t
    
    Vcmax25 = 60;%; % C3CRO;  BACOUR 2019 JGR
    
    Jmax25 = (2.59-0.035.*Temperature_dataset0).*Vcmax25; %Kattge and Knorr (2007);
    
    a_LL = 0.5;
    theta = 0.7; %von Caemmerer 2000
    Jmax = Jmax25.*exp(-((Temperature_dataset0-25)./18).^2); % Yin 2009； June 2004 Functional plant biology
    ETR_J = (a_LL.*alpha.*PAR_dataset0+Jmax-((a_LL.*alpha.*PAR_dataset0+Jmax).^2-4.*theta.*Jmax.*a_LL.*alpha.*PAR_dataset0).^(0.5))./(2.*theta);
        
    vwc_max = 0.40;
    vwc_min = 0.0875;
    fw = 10./3.*((vwc - vwc_min)./(vwc_max - vwc_min));
    fw (fw>1)=1;
    fw (fw<0)=0;
    % iteration
    Cc=0.7.*Ca;
    delta_Cc =70;
    i=0;
    
    while delta_Cc >0.1
        
        % calculate Tau_star  is the CO2 compensation point
        R = 8.314;% J mol K
        O = 230000; % Oxygen partial pressure in ubar (PPM)
        Sc_o25=2800;
        ESc_o = -24460;
        Sc_o =  Sc_o25.*exp((Temperature_dataset0-25).*ESc_o./(298.*R.*(Temperature_dataset0+273)));
        Tau_star = 0.5.*O./Sc_o;
        
        % calculate Vcmax and parameters in Ac
        Evcmax = 65330;  % Bernacchi  (2001) Plant cell environment
        R = 8.314;% J mol K
        Vcmax = Vcmax25.*exp((Temperature_dataset0-25).*Evcmax./(298.*R.*(Temperature_dataset0+273)));
        
        KmC25 = 270; % Yin 2009
        EKmC25=80990; % Yin 2009
        KmC = KmC25.*exp((Temperature_dataset0-25).*EKmC25./(298.*R.*(Temperature_dataset0+273)));
        
        KmO25 = 165000; % Yin 2009
        EKmO25 = 23720; % Yin 2009
        KmO = KmO25.*exp((Temperature_dataset0-25).*EKmO25./(298.*R.*(Temperature_dataset0+273)));
        
        Rd25 = 0.015.*Vcmax25; % Yin 2009
        ERd = 46390; % Yin 2009
        Rd = Rd25.*exp((Temperature_dataset0-25).*ERd./(298.*R.*(Temperature_dataset0+273)));
        
        Ac =(Cc-Tau_star).*Vcmax./(Cc+KmC.*(1+O./KmO)) -Rd; % Yin 2009  Eq2
        
        Aj =   ETR_J.*((Cc-Tau_star)./(4.*Cc+8.*Tau_star))-Rd;
        Anet = min([Aj Ac], [], 2);
        
        % Gc
        Go=0.01;
        a_1=11.0;
        Do=1.5;
        Cs = a_1./(a_1-1).*Cc;
        Gs = Go+(a_1.*fw.*Anet)./(Cs.*(1+VPD./Do));
        Gc = 0.64.*Gs;
        Cc_temp=Ca-Anet./Gc;
        % difference between Cc_temp and Cc
        delta_Cc = mean(abs(Cc-Cc_temp),'omitnan');
        Cc = Cc_temp;
        clear Cc_temp
        i=i+1;
    end
    
    % calculate Tau_star  is the CO2 compensation point
    R = 8.314;% J mol K
    O = 230000; % Oxygen partial pressure in ubar (PPM)
    Sc_o25=2800;
    ESc_o = -24460;
    Sc_o =  Sc_o25.*exp((Temperature_dataset0-25).*ESc_o./(298.*R.*(Temperature_dataset0+273)));
    Tau_star = 0.5.*O./Sc_o;
    
    % calculate Vcmax and parameters in Ac
    Evcmax = 65330;  % Bernacchi  (2001) Plant cell environment
    R = 8.314;% J mol K
    Vcmax = Vcmax25.*exp((Temperature_dataset0-25).*Evcmax./(298.*R.*(Temperature_dataset0+273)));
    
    KmC25 = 270; % Yin 2009
    EKmC25=80990; % Yin 2009
    KmC = KmC25.*exp((Temperature_dataset0-25).*EKmC25./(298.*R.*(Temperature_dataset0+273)));
    
    KmO25 = 165000; % Yin 2009
    EKmO25 = 23720; % Yin 2009
    KmO = KmO25.*exp((Temperature_dataset0-25).*EKmO25./(298.*R.*(Temperature_dataset0+273)));
    
    Rd25 = 0.015.*Vcmax25; % Yin 2009
    ERd = 46390; % Yin 2009
    Rd = Rd25.*exp((Temperature_dataset0-25).*ERd./(298.*R.*(Temperature_dataset0+273)));
    
    Ac =(Cc-Tau_star).*Vcmax./(Cc+KmC.*(1+O./KmO)) -Rd; % Yin 2009  Eq2
    
    ETRc_psii = (Ac+Rd)./((Cc-Tau_star)./(4.*Cc+8.*Tau_star)); % Yin 2009  Eq3a
    
    % find min between two ETRs
    ETR = min([ETR_J ETRc_psii], [], 2);
    YII = ETR./(PAR_dataset0.*alpha.*0.5);
    % normailzed between 0.1 and 0.8 to match with Tol 2014
    PSIIP0_max = 0.8;
    PSIIP0_min = 0.1;
    YII =  (PSIIP0_max-PSIIP0_min)*(YII-min(YII))/(max(YII)-min(YII)) + PSIIP0_min;
    x = 1-YII/PSIIP0_max; %  van der Tol 2014
    
    %     BACOUR 2019 Eq6
    a =16.042;
    b=5.74;
    c=2.167;
    d = -0.014;
    e=-0.00437;
    f=0.00576;
    Temperature = Temperature_dataset0;
    KNPQ = a.*x.^c.*(1+b)./(b+x.^c).*(exp(d.*Temperature_dataset0+e)./(PAR_dataset0.^f));
    KNPQ_x = a.*x.^c.*(1+b)./(b+x.^c);
    KNPQ_t = exp(d.*Temperature_dataset0+e);
    KNPQ_PAR = PAR_dataset0.^f;
    
    [m, ~] =find(KNPQ==0);
    KNPQ(m)=KNPQ(m-1);
    
    %     BACOUR 2019 Eq3
    KF=0.1;
    KD=0.9;
    phiF0_psii=0.02;
    phiF_psii =(KF./(KF+KD+KNPQ)).*(1-YII);    
    epsilon = phiF_psii./phiF0_psii;    
    NPQ = KNPQ;    
    Kp = YII./(1-YII).*(1+KNPQ);
    
    % phiF, phiP, phiD,phiN
    phi_N = KNPQ./(1+Kp+KNPQ);
    phi_F =  KF./(1+Kp+KNPQ);
    phi_D =  KD./(1+Kp+KNPQ);
    phi_P =  Kp./(1+Kp+KNPQ);
    
    %     BACOUR 2019 Eq1
    ChlF_PSII_ratio =0.00917.*epsilon./(0.00917.*epsilon+0.00561);
    fLP = 0.88;
    SIF_dataset0.SVD_O2A_PSII = ChlF_PSII_ratio.*SIF_dataset0.SVD_O2A;
    SIF_tot_O2A = SIF_dataset0.SVD_O2A_PSII./(fLP.*fesc); % Xinjie Liu AFM 2021
    
    % convert to full spectrum in ChlF at PS level
    [~,S,V] = svd(RC_PSII);
    SIF_PS_PSII_O2A = V(:,1)'.*(SIF_tot_O2A./V(760-640+1,1));
    SIF_PS_PSII_O2A_PPFD = PPFD_ChlF(640:850, SIF_PS_PSII_O2A);
    % GPP_SIF
    GPP = 0.25.*(Cc-Tau_star)./(Cc+2*Tau_star).*YII./(1-YII).*(1+NPQ).*10.*SIF_PS_PSII_O2A_PPFD;
    
    J_SIF  = YII./(1-YII).*(1+NPQ).*10.*SIF_PS_PSII_O2A_PPFD;
    GPP_LUE  =0.25.*(Cc-Tau_star)./(Cc+2*Tau_star).*phi_P.*PAR_dataset0.*SIF_dataset0.fPAR.*0.5 ;
    
    
    fPSII= ChlF_PSII_ratio;%SIF_dataset0.SVD_O2A_PSII./SIF_dataset0.SVD_O2A;
    
    SIF_dataset0 = [SIF_dataset0 table(fesc) table(SIF_tot_O2A) table(GPP) table(GPP_dataset0) table(PAR_dataset0)...
        table(Temperature_dataset0) table(YII) table(NPQ) table(fPSII) table(SIF_PS_PSII_O2A_PPFD)...
        table(phi_N) table(phi_F) table(phi_D) table(phi_P) table(Kp) table(KNPQ) table(epsilon)  table(J_SIF) table(GPP_LUE)];
    
    index2  = find( SIF_dataset0.GPP <=0 | SIF_dataset0.GPP_dataset0 <=0);
    SIF_dataset0(index2,:)=[];
    
    clear index2
    eval(['SIF_GPP_dataset_temp',num2str(ii),'=SIF_dataset0;'])
    
end


% organization
[C,ia,ic] = unique([SIF_GPP_dataset_temp1.Datetime
    SIF_GPP_dataset_temp2.Datetime
    SIF_GPP_dataset_temp3.Datetime
    SIF_GPP_dataset_temp4.Datetime
    SIF_GPP_dataset_temp5.Datetime
    SIF_GPP_dataset_temp6.Datetime]);

sz = [size(C,1),size(SIF_GPP_dataset_temp1,2)];
varNames = SIF_GPP_dataset_temp1.Properties.VariableNames;
varTypes = {'datetime','double','double','double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double'};
SIF_GPP_dataset1 = table('size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
SIF_GPP_dataset1.Datetime = C;
SIF_GPP_dataset1( ismember(C,SIF_GPP_dataset_temp1.Datetime),2:end) = SIF_GPP_dataset_temp1(:,2:end);
SIF_GPP_dataset1 = standardizeMissing(SIF_GPP_dataset1,0);

SIF_GPP_dataset2 = table('size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
SIF_GPP_dataset2.Datetime = C;
SIF_GPP_dataset2( ismember(C,SIF_GPP_dataset_temp2.Datetime),2:end) = SIF_GPP_dataset_temp2(:,2:end);
SIF_GPP_dataset2 = standardizeMissing(SIF_GPP_dataset2,0);

SIF_GPP_dataset3 = table('size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
SIF_GPP_dataset3.Datetime = C;
SIF_GPP_dataset3( ismember(C,SIF_GPP_dataset_temp3.Datetime),2:end) = SIF_GPP_dataset_temp3(:,2:end);
SIF_GPP_dataset3 = standardizeMissing(SIF_GPP_dataset3,0);

SIF_GPP_dataset4 = table('size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
SIF_GPP_dataset4.Datetime = C;
SIF_GPP_dataset4( ismember(C,SIF_GPP_dataset_temp4.Datetime),2:end) = SIF_GPP_dataset_temp4(:,2:end);
SIF_GPP_dataset4 = standardizeMissing(SIF_GPP_dataset4,0);

SIF_GPP_dataset5 = table('size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
SIF_GPP_dataset5.Datetime = C;
SIF_GPP_dataset5( ismember(C,SIF_GPP_dataset_temp5.Datetime),2:end) = SIF_GPP_dataset_temp5(:,2:end);
SIF_GPP_dataset5 = standardizeMissing(SIF_GPP_dataset5,0);

SIF_GPP_dataset6 = table('size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
SIF_GPP_dataset6.Datetime = C;
SIF_GPP_dataset6(ismember(C,SIF_GPP_dataset_temp6.Datetime),2:end) = SIF_GPP_dataset_temp6(:,2:end);
SIF_GPP_dataset6 = standardizeMissing(SIF_GPP_dataset6,0);

%% Table 2
run ('Table2_cal.m')

%% Figure4a DOY SIF760toc
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;
y_input =10.*mean([SIF_GPP_dataset1.SVD_O2A SIF_GPP_dataset2.SVD_O2A SIF_GPP_dataset3.SVD_O2A...
    SIF_GPP_dataset4.SVD_O2A SIF_GPP_dataset5.SVD_O2A SIF_GPP_dataset6.SVD_O2A ],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
ylim([0 2.6])
xlabel(['Month/year'],'Color','k')
ylabel(['SIF_T_O_C_\__7_6_0  (mW m^-^2 nm^-^1 sr^-^1)'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
%print('D:\Dropbox\Fluorescence_Liu_10\figure\temp\Figure3_YL_DOY_SIF760','-dtiff','-r300')
print('Figure4a','-djpeg','-r300')

%% Figure4b DOY GPPec
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.GPP_dataset0 SIF_GPP_dataset2.GPP_dataset0 SIF_GPP_dataset3.GPP_dataset0...
    SIF_GPP_dataset4.GPP_dataset0 SIF_GPP_dataset5.GPP_dataset0 SIF_GPP_dataset6.GPP_dataset0 ],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
ylim([0 max(y_input)])
xlabel(['Month/year'],'Color','k')
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('Figure4b','-djpeg','-r300')
%% Figure4c DOY NIRv
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.NDVI.*SIF_GPP_dataset1.NIR_v_O2A	 SIF_GPP_dataset2.NDVI.*SIF_GPP_dataset2.NIR_v_O2A SIF_GPP_dataset3.NDVI.*SIF_GPP_dataset3.NIR_v_O2A...
    SIF_GPP_dataset4.NDVI.*SIF_GPP_dataset4.NIR_v_O2A	 SIF_GPP_dataset5.NDVI.*SIF_GPP_dataset5.NIR_v_O2A SIF_GPP_dataset6.NDVI.*SIF_GPP_dataset6.NIR_v_O2A],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
ylim([0 1])
xlabel(['Month/year'],'Color','k')
ylabel(['NIR_v'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('Figure4c','-djpeg','-r300')

%% Figure4d DOY fAPAR
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.13]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));

x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.fPAR SIF_GPP_dataset2.fPAR SIF_GPP_dataset3.fPAR...
    SIF_GPP_dataset4.fPAR SIF_GPP_dataset5.fPAR SIF_GPP_dataset6.fPAR],2,'omitnan');
Hour_group_dataset0 = [SIF_GPP_dataset1 table(y_input)];
daily_temp_x = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','Datetime');
x_input2 = daily_temp_x.mean_Datetime;
daily_temp_y = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','y_input');
y_input2 = daily_temp_y.mean_y_input;
g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')

xtickformat('MM/yyyy')
xlabel(['Month/year'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
%set(gca,'Xticklabel',[])
ylabel(['{\itf}_A_P_A_R '])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'a'};
% text(datetime(max(x_input)+days(1)),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',32)
print('Figure4d','-djpeg','-r300')

%% Figure4e DOY Rred+Rnir
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.13]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.NIR_v_O2B SIF_GPP_dataset2.NIR_v_O2B SIF_GPP_dataset3.NIR_v_O2B...
    SIF_GPP_dataset4.NIR_v_O2B SIF_GPP_dataset5.NIR_v_O2B SIF_GPP_dataset6.NIR_v_O2B],2,'omitnan');
Hour_group_dataset0 = [SIF_GPP_dataset1 table(y_input)];
daily_temp_x = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','Datetime');
x_input2 = daily_temp_x.mean_Datetime;
daily_temp_y = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','y_input');
y_input2 = daily_temp_y.mean_y_input;
yyaxis  left
g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[255 177 177],...
    'MarkerFaceColor',1/255.*[255 177 177]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
ylim([0 0.82])
set(gca,'YTick',[0: 0.2 :1],'ycolor','r')
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
xlabel(['Month/year'],'Color','k')
ylabel(['{\itR}_6_8_0'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
hold off

yyaxis right

x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.NIR_v_O2A SIF_GPP_dataset2.NIR_v_O2A SIF_GPP_dataset3.NIR_v_O2A...
    SIF_GPP_dataset4.NIR_v_O2A SIF_GPP_dataset5.NIR_v_O2A SIF_GPP_dataset6.NIR_v_O2A],2,'omitnan');
Hour_group_dataset0 = [SIF_GPP_dataset1 table(y_input)];
daily_temp_x = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','Datetime');
x_input2 = daily_temp_x.mean_Datetime;
daily_temp_y = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','y_input');

y_input2 = daily_temp_y.mean_y_input;

g3 = plot(x_input,y_input,'*',...
    'MarkerSize',4,...
    'MarkerEdgeColor',1/255.*[255 177 255],...
    'MarkerFaceColor',1/255.*[255 177 255]);
hold on
g4 =plot(x_input2,y_input2,'s',...
    'MarkerSize',7,...
    'MarkerEdgeColor','m',...
    'MarkerFaceColor','m');
hold off
ylim([0 0.82])
set(gca,'YTick',[0: 0.2 :1],'ycolor','m')
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')

xtickformat('MM/yyyy')
xlabel(['Month/year'],'Color','k')
ylabel(['{\itR}_7_5_5'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on

legend([g1,g2,g3,g4], '{\itR}_6_8_0 Half hourly','{\itR}_6_8_0 Daily','{\itR}_7_5_5 Half hourly','{\itR}_7_5_5 Daily','Location','north','NumColumns',2,'FontSize',20)
legend('boxoff')
% txt = {'a'};
% text(datetime(max(x_input)+days(1)),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',32)
print('Figure4e','-djpeg','-r300')

%% Figure4f DOY NDVI
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.NDVI SIF_GPP_dataset5.NDVI],2,'omitnan');
%mean([SIF_GPP_dataset1.NDVI SIF_GPP_dataset2.NDVI SIF_GPP_dataset3.NDVI...
% SIF_GPP_dataset4.NDVI SIF_GPP_dataset5.NDVI SIF_GPP_dataset6.NDVI],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
%ylim([0 1500])
xlabel(['Month/year'],'Color','k')
ylabel(['NDVI'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('Figure4f','-djpeg','-r300')

%% Figure4g DOY Tair
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.Temperature_dataset0 SIF_GPP_dataset2.Temperature_dataset0 SIF_GPP_dataset3.Temperature_dataset0...
    SIF_GPP_dataset4.Temperature_dataset0 SIF_GPP_dataset5.Temperature_dataset0 SIF_GPP_dataset6.Temperature_dataset0 ],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
%ylim([0 1500])
xlabel(['Month/year'],'Color','k')
ylabel(['Air temperature (^oC)'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('Figure4g','-djpeg','-r300')

%% Figure4h DOY PAR
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.PAR_dataset0 SIF_GPP_dataset2.PAR_dataset0 SIF_GPP_dataset3.PAR_dataset0...
    SIF_GPP_dataset4.PAR_dataset0 SIF_GPP_dataset5.PAR_dataset0 SIF_GPP_dataset6.PAR_dataset0 ],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
% ylim([0 2000])
xlabel(['Month/year'],'Color','k')
ylabel(['PAR (\mumol m^-^2 s^-^1)'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('Figure4h','-djpeg','-r300')

%% Figure5 GPPsif GPPec half hourly+daily
ax1= figure('visible','on');
set(gcf,'Position',[50 50 1500 600]);
h1 = subplot(1,2,1);
set(gca,'position', [0.08 0.18 0.40 0.8]);   % [x0 y0 width height]

x_input_temp = mean([SIF_GPP_dataset1.GPP SIF_GPP_dataset2.GPP SIF_GPP_dataset3.GPP...
    SIF_GPP_dataset4.GPP SIF_GPP_dataset5.GPP SIF_GPP_dataset6.GPP ],2,'omitnan');
y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0 SIF_GPP_dataset2.GPP_dataset0 SIF_GPP_dataset3.GPP_dataset0...
    SIF_GPP_dataset4.GPP_dataset0 SIF_GPP_dataset5.GPP_dataset0 SIF_GPP_dataset6.GPP_dataset0 ],2,'omitnan');

t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input = SIF_GPP_dataset1.PAR_dataset0;
x_input =t1;
y_input = t2;

y=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
x=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
h(1) = scatter(x,y,50,color_bar,'o','filled');
hold on

h1= colorbar('eastoutside');
colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
h1.Limits = [0 	1900];
h1.Ticks = [0:300:1800];

xlim([0 75])
set(gca,'XTick',[0: 20 :60],'ycolor','k')
ylim([0 75])
set(gca,'YTick',[0: 20 :60],'ycolor','k')
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
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);

box on
% txt2 = { 'Half hourly'} ;
% text(0.4*max(xlim),1.08*max(ylim),txt2,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
txt = {'a'};
text(0.92*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',40)



h2 = subplot(1,2,2);
set(gca,'position', [0.62 0.18 0.36 0.8]);   % [x0 y0 width height]

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input_temp)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input_temp');
x_input_daily = GPP_daily_temp.mean_x_input_temp;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input_temp)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input_temp');
y_input_daily = GPP_dataset0_daily_temp.mean_y_input_temp;

x_input = x_input_daily;
y_input = y_input_daily;

g2 = plot(x_input,y_input,'k^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
%set(get(get(g2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
if max(x_input) > max(y_input)
    xlim([0 1.1*max(x_input)])
    ylim([0 1.1*max(x_input)])
else
    xlim([0 1.1*max(y_input) ])
    ylim([0 1.1*max(y_input) ])
end

xlim([0 45])
set(gca,'XTick',[0: 10 :40],'ycolor','k')
xlabel(['GPP_S_I_F (\mumol m^-^2 s^-^1)'],'Color','k')

ylim([0 45])
set(gca,'YTick',[0: 10 :40],'ycolor','k')
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
hold on
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
plot(xrange, xrange,'r-.','LineWidth',2) % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
plot(xrange,ypred,'k-','LineWidth',1.5)
boundedline(xrange, ypred, delta,'alpha');
% plot(xrange,delta_t,'r--','LineWidth',1)
% plot(xrange,ypred,'k-','LineWidth',1.5)
hold off
lgd = legend('Daily','1:1 line','Linear fit','Location','northwest','NumColumns',1,'FontSize',22);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(RMSEtest(x_input,y_input),-4);
rrmse = 100.*RMSE2./(max(x_input)-min(x_input));
title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
%set(gca,'LooseInset',get(gca,'TightInset'));
box on
% txt2 = { 'Half hourly'} ;
% text(0.4*max(xlim),1.08*max(ylim),txt2,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
txt = {'b'};
text(0.92*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',40)
print('Figure5','-djpeg','-r300')



%% Figure6a DOY fPSII
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
%ha = tight_subplot(2,2,[0.08 0.09],[.09 .02],[0.1 0.05]);
%axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.fPSII SIF_GPP_dataset2.fPSII SIF_GPP_dataset3.fPSII...
    SIF_GPP_dataset4.fPSII SIF_GPP_dataset5.fPSII SIF_GPP_dataset6.fPSII],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
ylim([0.3 0.8])
set(gca,'YTick',[0.3: 0.1 :0.75],'ycolor','k')
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
xlabel(['Month/year'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
%set(gca,'Xticklabel',[])
ylabel(['{\itf}_P_S_I_I'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'c'};
% text(datetime(max(x_input)+days(5)),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',32)
print('Figure6a','-djpeg','-r300')

%% Figure6b Tair fPSII
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
%ha = tight_subplot(2,2,[0.08 0.09],[.09 .02],[0.1 0.05]);
%axes(ha(1));
x_input = SIF_GPP_dataset1.Temperature_dataset0;
y_input =mean([SIF_GPP_dataset1.fPSII SIF_GPP_dataset2.fPSII SIF_GPP_dataset3.fPSII...
    SIF_GPP_dataset4.fPSII SIF_GPP_dataset5.fPSII SIF_GPP_dataset6.fPSII],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
ylim([0.3 0.8])

xlabel(['Air temperature (^oC)'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
%set(gca,'Xticklabel',[])
ylabel(['{\itf}_P_S_I_I'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'c'};
% text(datetime(max(x_input)+days(5)),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',32)
print('Figure6b','-djpeg','-r300')

%% Figure6c PAR fPSII
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
%ha = tight_subplot(2,2,[0.08 0.09],[.09 .02],[0.1 0.05]);
%axes(ha(1));
x_input = SIF_GPP_dataset1.PAR_dataset0;
y_input =mean([SIF_GPP_dataset1.fPSII SIF_GPP_dataset2.fPSII SIF_GPP_dataset3.fPSII...
    SIF_GPP_dataset4.fPSII SIF_GPP_dataset5.fPSII SIF_GPP_dataset6.fPSII],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
ylim([0.3 0.8])
xlabel(['PAR (\mumol m^-^2 s^-^1)'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
%set(gca,'Xticklabel',[])
ylabel(['{\itf}_P_S_I_I'])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'c'};
% text(datetime(max(x_input)+days(5)),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',32)
print('Figure6c','-djpeg','-r300')

%% Figure6d DOY fesc_P-C
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

% ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
% set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;
y_input =SIF_GPP_dataset1.fesc;%mean([SIF_GPP_dataset1.fesc SIF_GPP_dataset2.fesc SIF_GPP_dataset3.fesc...
% SIF_GPP_dataset4.fesc SIF_GPP_dataset5.fesc SIF_GPP_dataset6.fesc],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
ylim([0 1])
xlabel(['Month/year'],'Color','k')
ylabel(['{\itf}_e_s_c_\__P_-_C'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('Figure6d','-djpeg','-r300')


%% Figure7 phi_4_daily
GPP_daily_temp = groupsummary(SIF_GPP_dataset1,'Datetime','dayofyear','mean','phi_F');
phi_F_daily = GPP_daily_temp.mean_phi_F;

GPP_daily_temp = groupsummary(SIF_GPP_dataset1,'Datetime','dayofyear','mean','phi_D');
phi_D_daily = GPP_daily_temp.mean_phi_D;

GPP_daily_temp = groupsummary(SIF_GPP_dataset1,'Datetime','dayofyear','mean','phi_N');
phi_N_daily = GPP_daily_temp.mean_phi_N;

GPP_daily_temp = groupsummary(SIF_GPP_dataset1,'Datetime','dayofyear','mean','phi_P');
phi_P_daily = GPP_daily_temp.mean_phi_P;

GPP_daily_temp = groupsummary(SIF_GPP_dataset1,'Datetime','dayofyear','mean','Datetime');
Datetime_daily = GPP_daily_temp.mean_Datetime;

phi_4_daily = [phi_F_daily phi_D_daily phi_N_daily phi_P_daily];


x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.PAR_dataset0 SIF_GPP_dataset2.PAR_dataset0 SIF_GPP_dataset3.PAR_dataset0...
    SIF_GPP_dataset4.PAR_dataset0 SIF_GPP_dataset5.PAR_dataset0 SIF_GPP_dataset6.PAR_dataset0 ],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
PAR = GPP_dataset0_daily_temp.mean_y_input;



x_input = SIF_GPP_dataset1.Datetime;
y_input =mean([SIF_GPP_dataset1.Temperature_dataset0 SIF_GPP_dataset2.Temperature_dataset0 SIF_GPP_dataset3.Temperature_dataset0...
    SIF_GPP_dataset4.Temperature_dataset0 SIF_GPP_dataset5.Temperature_dataset0 SIF_GPP_dataset6.Temperature_dataset0 ],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
Temperature = GPP_dataset0_daily_temp.mean_y_input;


% figure 6_2
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.21 .025],[0.11 0.03]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 13]);

h = area(x_input2,phi_4_daily,'LineStyle','-');
h(1).FaceColor = 'r';
h(2).FaceColor = '[0.4 0.4 0.4]';
h(3).FaceColor = 'b';
h(4).FaceColor = '[0 0.5 0]';

ylim([0 1.2])
set(gca,'YTick',[0: 0.2 :1.1],'ycolor','k')

x_input = Datetime_daily;
xlim([min(x_input)-days(1) max(x_input)+days(1)])
%set(gca,'XTick',[min(x_input)-days(30): calmonths(1): max(x_input)],'xcolor','k')

xtickformat('MM/yyyy')
xlabel(['Month/year'],'Color','k')
ylabel(['Yields \Phi (Daily)'])
legend('\Phi_F','\Phi_D','\Phi_N','\Phi_P','Location','north','NumColumns',4,'FontSize',16')
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',16);
set(gca,'LooseInset',get(gca,'TightInset'));
print('Figure7_1','-djpeg','-r300')

figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.06 .26],[0.03 0.2]);
set(gcf,'Units','centimeters', 'Position',[5 5 14 10]);
yyaxis right
h = plot(x_input2,phi_F_daily,'r-','LineWidth',2);
ylim([0 0.03])
set(gca,'YTick',[0.0: 0.01 :0.03],'ycolor','k')

x_input = Datetime_daily;
xlim([min(x_input)-days(1) max(x_input)+days(1)])
xtickformat('MM/yyyy')
%xlabel(['Month/year'],'Color','k');
set(gca,'XTick',[],'xcolor','k')
%ylabel(['\Phi_F ( Daily)'])

set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',40);
set(gca,'LooseInset',get(gca,'TightInset'));
title('\Phi_F (Daily)')
print('Figure7_2','-djpeg','-r300')

%% FigureS3 fesc_NDVI
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);
axes(ha(1));
%ha = tight_subplot(2,2,[0.08 0.09],[.09 .02],[0.1 0.05]);
%axes(ha(1));
x_input = SIF_GPP_dataset1.NDVI;
y_input =SIF_GPP_dataset1.fesc;

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
ylim([0 1])

xlabel(['NDVI'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
%set(gca,'Xticklabel',[])
ylabel(['{\itf}_e_s_c_\__P_-_C'],'Color','k')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'c'};

% text(datetime(max(x_input)+days(5)),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',32)
print('FigureS3','-djpeg','-r300')

%% FigureS4a DOY-Kn
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

% ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
% set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;

y_input =mean([SIF_GPP_dataset1.KNPQ SIF_GPP_dataset2.KNPQ SIF_GPP_dataset3.KNPQ...
    SIF_GPP_dataset4.KNPQ SIF_GPP_dataset5.KNPQ SIF_GPP_dataset6.KNPQ],2,'omitnan');

% y_input =mean([SIF_GPP_dataset1.phi_F SIF_GPP_dataset2.phi_F SIF_GPP_dataset3.phi_F...
%     SIF_GPP_dataset4.phi_F SIF_GPP_dataset5.phi_F SIF_GPP_dataset6.phi_F],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
% ylim([0 1])
xlabel(['Month/year'],'Color','k')
ylabel(['{\itK}_N'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('FigureS4a','-djpeg','-r300')
%% FigureS4b DOY-Kp
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

axes(ha(1));
x_input = SIF_GPP_dataset1.Datetime;

y_input =mean([SIF_GPP_dataset1.Kp SIF_GPP_dataset2.Kp SIF_GPP_dataset3.Kp...
    SIF_GPP_dataset4.Kp SIF_GPP_dataset5.Kp SIF_GPP_dataset6.Kp],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
xlim([min(x_input)-days(5) max(x_input)+days(5)])
set(gca,'XTick',[min(x_input)+days(17): calmonths(2): max(x_input)],'xcolor','k')
xtickformat('MM/yyyy')
% ylim([0 1])
xlabel(['Month/year'],'Color','k')
ylabel(['{\itK}_P'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('FigureS4b','-djpeg','-r300')
%% FigureS4c Tair-Kn
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

axes(ha(1));
x_input = SIF_GPP_dataset1.Temperature_dataset0;

y_input =mean([SIF_GPP_dataset1.KNPQ SIF_GPP_dataset2.KNPQ SIF_GPP_dataset3.KNPQ...
    SIF_GPP_dataset4.KNPQ SIF_GPP_dataset5.KNPQ SIF_GPP_dataset6.KNPQ],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
% set(gca,'XTick',[min(x_input): calmonths(2): max(x_input)+calmonths(1) ],'xcolor','k')
% xlim([min(x_input)-days(10) max(x_input)+days(15)])
% xtickformat('MM/yyyy')
% % ylim([0 1])
xlabel(['Air temperature (^oC)'],'Color','k')
ylabel(['{\itK}_N'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('FigureS4c','-djpeg','-r300')
%% FigureS4d Tair-Kn
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

axes(ha(1));
x_input = SIF_GPP_dataset1.Temperature_dataset0;

y_input =mean([SIF_GPP_dataset1.Kp SIF_GPP_dataset2.Kp SIF_GPP_dataset3.Kp...
    SIF_GPP_dataset4.Kp SIF_GPP_dataset5.Kp SIF_GPP_dataset6.Kp],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
% set(gca,'XTick',[min(x_input): calmonths(2): max(x_input)+calmonths(1) ],'xcolor','k')
% xlim([min(x_input)-days(10) max(x_input)+days(15)])
% xtickformat('MM/yyyy')
% % ylim([0 1])
xlabel(['Air temperature (^oC)'],'Color','k')
ylabel(['{\itK}_P'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('FigureS4d','-djpeg','-r300')

%% FigureS4e PAR-Kn
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

axes(ha(1));
x_input = SIF_GPP_dataset1.PAR_dataset0;

y_input =mean([SIF_GPP_dataset1.KNPQ SIF_GPP_dataset2.KNPQ SIF_GPP_dataset3.KNPQ...
    SIF_GPP_dataset4.KNPQ SIF_GPP_dataset5.KNPQ SIF_GPP_dataset6.KNPQ],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
% set(gca,'XTick',[min(x_input): calmonths(2): max(x_input)+calmonths(1) ],'xcolor','k')
% xlim([min(x_input)-days(10) max(x_input)+days(15)])
% xtickformat('MM/yyyy')
% % ylim([0 1])
xlabel(['PAR (\mumol m^-^2 s^-^1)'],'Color','k')
ylabel(['{\itK}_N'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('FigureS4e','-djpeg','-r300')

%% FigureS4f PAR-Kn
figure('visible','on');
ha = tight_subplot(1,1,[0.08 0.08],[.18 .02],[0.16 0.12]);
set(gcf,'Units','centimeters', 'Position',[5 5 22 16]);

axes(ha(1));
x_input = SIF_GPP_dataset1.PAR_dataset0;

y_input =mean([SIF_GPP_dataset1.Kp SIF_GPP_dataset2.Kp SIF_GPP_dataset3.Kp...
    SIF_GPP_dataset4.Kp SIF_GPP_dataset5.Kp SIF_GPP_dataset6.Kp],2,'omitnan');

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input');
x_input2 = GPP_daily_temp.mean_x_input;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input');
y_input2 = GPP_dataset0_daily_temp.mean_y_input;

g1 = plot(x_input,y_input,'.',...
    'MarkerSize',12,...
    'MarkerEdgeColor',1/255.*[167 211 255],...
    'MarkerFaceColor',1/255.*[167 211 255]);
hold on
g2 =plot(x_input2,y_input2,'^',...
    'MarkerSize',7,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
hold off
% set(gca,'XTick',[min(x_input): calmonths(2): max(x_input)+calmonths(1) ],'xcolor','k')
% xlim([min(x_input)-days(10) max(x_input)+days(15)])
% xtickformat('MM/yyyy')
% % ylim([0 1])
xlabel(['PAR (\mumol m^-^2 s^-^1)'],'Color','k')
ylabel(['{\itK}_P'],'Color','k')
lgd = legend('Half hourly','Daily','Location','northwest','NumColumns',1,'FontSize',24);
legend('boxoff')
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'(a)'};
% text(datetime('2021/06/10'),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
print('FigureS4f','-djpeg','-r300')

%% FigureS5 SIF GPPec half hourly
ax1= figure('visible','on');
set(gcf,'Position',[50 50 1500 600]);
h1 = subplot(1,2,1);
set(gca,'position', [0.08 0.18 0.40 0.8]);   % [x0 y0 width height]


x_input_temp = 10.*mean([SIF_GPP_dataset1.SVD_O2A SIF_GPP_dataset2.SVD_O2A SIF_GPP_dataset3.SVD_O2A...
    SIF_GPP_dataset4.SVD_O2A SIF_GPP_dataset5.SVD_O2A SIF_GPP_dataset6.SVD_O2A ],2,'omitnan');
y_input_temp = mean([SIF_GPP_dataset1.GPP_dataset0 SIF_GPP_dataset2.GPP_dataset0 SIF_GPP_dataset3.GPP_dataset0...
    SIF_GPP_dataset4.GPP_dataset0 SIF_GPP_dataset5.GPP_dataset0 SIF_GPP_dataset6.GPP_dataset0 ],2,'omitnan');

%fitlm(x_input_temp,y_input_temp,'Intercept',false)

%sum(abs(y_input_temp./x_input_temp-1),'omitnan')/size(y_input_temp,1);


t1 = reshape(x_input_temp,[size(x_input_temp,1)*size(x_input_temp,2),1]);
t2 = reshape(y_input_temp,[size(y_input_temp,1)*size(y_input_temp,2),1]);

color_input = SIF_GPP_dataset1.PAR_dataset0;
x_input =t1;
y_input = t2;

y=reshape(y_input,size(y_input,1)*size(y_input,2),1);% may need to change the 1 or 2
x=repmat(x_input,size(y_input,2),1);
color_input_temp=repmat(color_input,size(y_input,2),1);
color_bar=reshape(color_input_temp,size(y_input,1)*size(y_input,2),1);
clear color_input_temp;
h(1) = scatter(x,y,50,color_bar,'d','filled');
hold on
%h(1) = plot(nan,nan,'kd','MarkerSize',20);
%set(get(get(g1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(h1, 'CLim', [min(color_input) max(color_input)]);
h1= colorbar('eastoutside');
colormap(jet);
h1.Label.String = 'PAR (\mumol m^-^2 s^-^1)';
h1.Limits = [0 	1900];
h1.Ticks = [0:300:1800];


xlim([0 2.5])
set(gca,'XTick',[0: 0.5 :2.5],'ycolor','k')
xlabel(['SIF_T_O_C_\__7_6_0  (mW m^-^2 nm^-^1 sr^-^1)'],'Color','k')

ylim([0 75])
set(gca,'YTick',[0: 20 :60],'ycolor','k')
%set(gca,'XTick',[]);
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');

model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
%plot(xrange, xrange,'r-.','LineWidth',2) % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
h(2) = plot(xrange,ypred,'k-','LineWidth',2);
boundedline(xrange, ypred, delta,'alpha');
% plot(xrange,delta_t,'r--','LineWidth',1)
% plot(xrange,ypred,'k-','LineWidth',1.5)
hold off
lgd = legend(h(1:2),'Half hourly','Linear fit','Location','northwest','NumColumns',1,'FontSize',22);
legend('boxoff')
R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(model_temp.RMSE ,-4);
rrmse = 100.*RMSE2./(max(model_temp.Fitted)-min(model_temp.Fitted));

%title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]});
title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});

if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
%set(gca,'LooseInset',get(gca,'TightInset'));
box on
% txt2 = { 'Half hourly'} ;
% text(0.4*max(xlim),1.08*max(ylim),txt2,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
txt = {'a'};
text(0.92*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',40)

h2 = subplot(1,2,2);
set(gca,'position', [0.62 0.18 0.36 0.8]);   % [x0 y0 width height]

Hour_group_dataset0 = [SIF_GPP_dataset1 table(x_input_temp)];
GPP_daily_temp = groupsummary(Hour_group_dataset0,'Datetime','dayofyear','mean','x_input_temp');
x_input_daily = GPP_daily_temp.mean_x_input_temp;

Hour_group_dataset1 = [SIF_GPP_dataset1 table(y_input_temp)];
GPP_dataset0_daily_temp = groupsummary(Hour_group_dataset1,'Datetime','dayofyear','mean','y_input_temp');
y_input_daily = GPP_dataset0_daily_temp.mean_y_input_temp;

x_input = x_input_daily;
y_input = y_input_daily;

g2 = plot(x_input,y_input,'ks',...
    'MarkerSize',7,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
xlim([0 1.4])
set(gca,'XTick',[0: 0.3 :1.5],'xcolor','k')
xlabel(['SIF_T_O_C_\__7_6_0  (mW m^-^2 nm^-^1 sr^-^1)'],'Color','k')

ylim([0 45])
set(gca,'YTick',[0: 10 :40],'ycolor','k')
ylabel(['GPP_E_C (\mumol m^-^2 s^-^1)'])
set(gca,'ycolor','k');
hold on
model_temp = fitlm(x_input,y_input);
SE = diag(sqrt(model_temp.CoefficientCovariance));
xrange = [min(x_input):(max(x_input)-min(x_input))./1000:max(x_input)]';
%plot(xrange, xrange,'r-.','LineWidth',2) % 1:1 line
[ypred,delta_t] = predict(model_temp,xrange,'Alpha',0.05,'Simultaneous',false);
delta = abs(delta_t-ypred);
plot(xrange,ypred,'k-','LineWidth',1.5)
boundedline(xrange, ypred, delta,'alpha');
% plot(xrange,delta_t,'r--','LineWidth',1)
% plot(xrange,ypred,'k-','LineWidth',1.5)
hold off
lgd = legend('Daily','Linear fit','Location','northwest','NumColumns',1,'FontSize',22);
legend('boxoff')

R22=roundn(model_temp.Rsquared.Ordinary ,-4);
RMSE2=roundn(model_temp.RMSE ,-4);
rrmse = 100.*RMSE2./(max(model_temp.Fitted)-min(model_temp.Fitted));

title(lgd,{['{\it{R}}^2 = ',num2str(R22,'%6.2f')]; ['RMSE = ',num2str(RMSE2,'%6.2f')]; ['rRMSE = ',num2str(rrmse,'%6.2f'),'%']});
if model_temp.Coefficients.pValue(1)<0.05
    txt = {...
        ['y-intercept = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
else
    txt = {...
        ['y-intercept^* = ',num2str(model_temp.Coefficients.Estimate(1),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(1),'%10.2g')],...
        ['slope = ',num2str(model_temp.Coefficients.Estimate(2),'%10.2g'),' ± ',num2str(model_temp.Coefficients.SE(2),'%10.2g')],...
        };
end
text(0.5*max(xlim),0.10*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
hold off
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
%set(gca,'LooseInset',get(gca,'TightInset'));
box on
% txt2 = { 'Half hourly'} ;
% text(0.4*max(xlim),1.08*max(ylim),txt2,'Fontname','times new Roman','FontWeight','bold','FontSize',18)
txt = {'b'};
text(0.92*max(xlim),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',40)

print('FigureS5','-djpeg','-r300')

close all

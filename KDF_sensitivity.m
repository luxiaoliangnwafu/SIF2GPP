clc;clear;
load RSE_dataset.mat
%% GPP calculation

R2_KDF = [];
RMSE_KDF= [];
RRMSE_KDF = [];
KDF_VALUES = [9:0.1:19]';
for jj = 1:size (KDF_VALUES,1)
    
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
        
        KDF=KDF_VALUES(jj);
        
        KF = 1/(KDF+1);
        KD = 1-KF;
        
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
        GPP = 0.25.*(Cc-Tau_star)./(Cc+2*Tau_star).*YII./(1-YII).*(1+NPQ).*(1+KDF).*SIF_PS_PSII_O2A_PPFD;
        
        J_SIF  = YII./(1-YII).*(1+NPQ).*(1+KDF).*SIF_PS_PSII_O2A_PPFD;
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
    
    model_temp = fitlm(x_input,y_input);
    
    R22=roundn(model_temp.Rsquared.Ordinary ,-4);
    RMSE2=roundn(RMSEtest(x_input,y_input),-4);
    rrmse = 100.*RMSE2./(max(x_input)-min(x_input));
    
    
    
    R2_KDF = [R2_KDF; R22];
    RMSE_KDF= [RMSE_KDF;RMSE2];
    RRMSE_KDF = [RRMSE_KDF; rrmse ];
    
    
    
end

%% Figure S6 KDF_sensitivity
figure('visible','on');
set(gcf,'Units','centimeters', 'Position',[4 5 22 16]);

x_input = KDF_VALUES;
y_input =R2_KDF;

x_input2 =KDF_VALUES;
y_input2 = RMSE_KDF;

yyaxis left
g1 = plot(x_input,y_input,'-','LineWidth',3,....
    'MarkerSize',5,...
    'MarkerEdgeColor','#0072BD',...
    'MarkerFaceColor','#0072BD');
xlim([9 19])
ylim([0.6 1.0])
ylabel(['{\it{R}}^2 GPP_E_C - GPP_S_I_F'])

set(gca,'XTick',[9: 2 :19],'xcolor','k')
set(gca,'YTick',[0.6: 0.1 :1])

yyaxis right
g2 =plot(x_input2,y_input2,'--','LineWidth',3,...
    'MarkerSize',5,...
    'MarkerEdgeColor','#D95319',...
    'MarkerFaceColor','#D95319');
xlim([9 19])
ylim([5 10])
xlabel(['K_D_F'],'Color','k')
ylabel(['RMSE  GPP_E_C - GPP_S_I_F'])

    lgd = legend('{\it{R}}^2','RMSE','Location','northwest','NumColumns',2,'FontSize',24);
legend('boxoff')
%set(gca,'Xticklabel',[])
set(gca,'Fontname','times new Roman','FontWeight','bold','FontSize',24);
set(gca,'LooseInset',get(gca,'TightInset'));
grid on
% txt = {'c'};
% text(datetime(max(x_input)+days(5)),0.95*max(ylim),txt,'Fontname','times new Roman','FontWeight','bold','FontSize',32)
print('FigureS6_KDF','-djpeg','-r500')











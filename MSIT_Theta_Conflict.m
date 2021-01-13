%% This script is to analyze theta evoked during high v low conflict trials of MSIT with no stimulation
clear all
MainPath='D:\BasuDDrive\MSIT paper\Neural Data\Theta I-C analysis\';
Subject='P09';
stimType='No_Stim';
load([MainPath,Subject,'PSD',stimType,'.mat']);

%% Pre screening for channel selection: 1) Finding channels that have enhanced theta after image onset in NoStim block trials, 2) Use all PFC channels
% 1. Use Theta pre and post image 

Id_freq=find(ft_freq.freq<9 & ft_freq.freq>3);
Id_time0=find(ft_freq.time>-0.6 & ft_freq.time<0);
Id_time1=find(ft_freq.time>0.2 & ft_freq.time<=1.2);

Theta0=nanmean(nanmean(ft_freq.powspctrm(:,:,Id_freq,Id_time0),4),3);
Theta=squeeze(nanmean(ft_freq.powspctrm(:,:,Id_freq,Id_time1),3));

H=zeros(size(Theta0,2),length(Id_time1));
for ch=1:size(Theta0,2)
    for tt=1:length(Id_time1)
        [~,~,ci]=ttest(log(Theta(:,ch,tt)),log(Theta0(:,ch)),'alpha',0.01);
        if(ci(1)>0 && ci(2)>0)
            H(ch,tt)=1;
        end
    end
end

Theta_channels=find(sum(H,2)>0);

% 2. Use all of PFC channels, specify the channel numbers
Pfc_channels=[1:41,51:65,99:136,150:161]';

%% Theta Analysis for High>Low conflict trials in the No-Stim blocks (NS1 trials)
Trials_C=find(~isnan(TrialDet(:,12)) & TrialDet(:,14)==1);
Trials_I=find(~isnan(TrialDet(:,12)) & TrialDet(:,14)==2);

% create new ft_data for C and I trials
ft_data_C=ft_data3;
ft_data_C.trial=ft_data3.trial(Trials_C);
ft_data_C.time=ft_data3.time(Trials_C);
ft_data_C.sampleinfo=ft_data3.sampleinfo(Trials_C,:);

ft_data_I=ft_data3;
ft_data_I.trial=ft_data3.trial(Trials_I);
ft_data_I.time=ft_data3.time(Trials_I);
ft_data_I.sampleinfo=ft_data3.sampleinfo(Trials_I,:);

% subtract ERP from I and C
addpath 'D:\C Drive\matlab codes\fieldtrip-master'

erp=ft_data_C.trial{1};
for ii=2:length(Trials_C)
erp=erp+ft_data_C.trial{ii};
end
ErpC=erp./length(Trials_C);

for cc=1:length(Trials_C)
   ft_data_C.trial{cc}= ft_data_C.trial{cc}-ErpC;
end

erp=ft_data_I.trial{1};
for ii=2:length(Trials_I)
erp=erp+ft_data_I.trial{ii};
end
ErpI=erp./length(Trials_I);

for ii=1:length(Trials_I)
   ft_data_I.trial{ii}= ft_data_I.trial{ii}-ErpI;
end

% Recalculate PSD 
cfg = [];
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.method = 'wavelet';
cfg.foi = [4:8];
cfg.toi = [-1.9:0.01:3.9];

ft_freqC = ft_freqanalysis(cfg,ft_data_C);
ft_freqI = ft_freqanalysis(cfg,ft_data_I);

Id_time0=find(ft_freqI.time>-0.6 & ft_freqI.time<0);
Id_time1=find(ft_freqI.time>0.1 & ft_freqI.time<=1.2);

% Use the PFC channels that have enhanced Theta post image onset
channels=intersect(Theta_channels,Pfc_channels);

Theta_C0=nanmean(nanmean(ft_freqC.powspctrm(:,channels,:,Id_time0),4),3);
Theta_C=nanmean(nanmean(ft_freqC.powspctrm(:,channels,:,Id_time1),4),3);

Theta_I0=nanmean(nanmean(ft_freqI.powspctrm(:,channels,:,Id_time0),4),3);
Theta_I=nanmean(nanmean(ft_freqI.powspctrm(:,channels,:,Id_time1),4),3);

save([loadPath,Subject,'PSD',stimType,'.mat'],'Trials_C','Trials_I','ThetaC','ThetaC0','ThetaI','ThetaI0','-append');

% Regression analysis
clear all
% Structure data
filename=dir('*.mat');
LogTheta=[];
Conflict=[];
Subject=[];
Region=[];
State=[];

for ff=1:length(filename)
    load(filename(ff).name,'theta_ictal','ThetaC','ThetaC0','ThetaI','ThetaI0','ParcellationValues','channels','score');
    ch=setdiff(1:size(ThetaC,2),theta_ictal);
    channels=channels(ch);
    Theta1=log(ThetaC(:,ch)./ThetaC0(:,ch));
    Theta2=log(ThetaI(:,ch)./ThetaI0(:,ch));
    
    sub=strings(size(Theta1,1)+size(Theta2,1),1);
    sub(:)=['P',num2str(ff)];
    for mi=1:length(channels)
        R=ParcellationValues(channels(mi),8);
        LogTheta=[LogTheta;Theta1(:,mi);Theta2(:,mi)];
        Conflict=[Conflict;zeros(length(Theta1(:,mi)),1);ones(length(Theta2(:,mi)),1)];
        Subject=[Subject;sub];
        Region=[Region;R.*ones(size(Theta1,1)+size(Theta2,1),1)];
        State=[State;sum(score).*ones(size(Theta1,1)+size(Theta2,1),1)];
   end
    
end

State_median=median(unique(State(~isnan(State))));
State_binary=zeros(size(State));
State_binary(State>State_median)=1;
State_binary(isnan(State))=nan;

% We are fitting GLME for each unique region
Region(Region==8)=6; % Both 6 and 8 encode vlPFC
R=unique(Region);
M=cell(length(R),1);
pval=[];

for rr=1:length(R)
    rid=find(Region==R(rr));
    L(rr)=length(rid);
    dataTable = table(LogTheta(rid),categorical(Conflict(rid)),categorical(State_binary(rid)),Subject(rid),...
            'VariableNames',{'LogTheta','Conflict','State','Subject'});
    M{rr}=fitglme(dataTable,'LogTheta ~ Conflict + State + (1|Subject)','Link','identity');
    pval(rr,1:2)=M{rr,1}.Coefficients.pValue(2:3);
end

[~,~,padj]=fdr(pval(:,1));

dataTable = table(LogTheta,categorical(Conflict),State,categorical(Region),Subject,...
            'VariableNames',{'LogTheta','Conflict','State','Region','Subject'});
Md=fitglme(dataTable,'LogTheta ~ Conflict + State + Region + (1|Subject)','Link','identity');

% Boxplot (Fig 2B)
Rid=[1,2,6,8,9]; % Plotting a subset of the regions
R=R(Rid);

x=[1,4,7,10,13];
for rr=1:length(R)
    rid=find(Region==R(rr));
    rtheta=LogTheta(rid);
    rconflict=Conflict(rid);
    bplot(rtheta(rconflict==0),x(rr),'std','width',0.6);hold on
    bplot(rtheta(rconflict==1),x(rr)+0.75,'std','width',0.6);
end
box off
xlim([-1 15])
ylim([-0.7 1.3])

xtlabel=[{'dlPFC'},{'dmPFC'},{'vlPFC'},{'dACC'},{'PPC'}]
set(gca,'xtick',x+0.4,'xticklabel',xtlabel)
legend('Non-conflict','Conflict')


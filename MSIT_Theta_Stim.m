%% This script is to analyze theta evoked during stimulation blocks of MSIT 
clear all
MainPath='E:\MSIT data for NatureBME\PSD\';
Subject='P11';
stimType='Stim1';
load([MainPath,'PSD_',Subject,'_',stimType,'.mat']);

%% Pre screening: Finding channels that have enhanced theta after image onset in NoStim block trials
% 1. Use Theta pre and post image 
Trials_NS0=find(strcmp(StimLocation,'None')==1 & ~isnan(TrialDet(:,12)));
Id_freq=find(ft_freq.freq<9 & ft_freq.freq>3);
Id_time0=find(ft_freq.time>-0.6 & ft_freq.time<0);
Id_time1=find(ft_freq.time>0.1 & ft_freq.time<=1.4);

Theta0=nanmean(nanmean(ft_freq.powspctrm(Trials_NS0,:,Id_freq,Id_time0),4),3);
Theta=squeeze(nanmean(ft_freq.powspctrm(Trials_NS0,:,Id_freq,Id_time1),3));

H=zeros(size(Theta0,2),length(Id_time1));
for ch=1:size(Theta0,2)
    for tt=1:length(Id_time1)
        [~,~,ci]=ttest(log(Theta(:,ch,tt)),log(Theta0(:,ch)),'alpha',0.05);
        if(ci(1)>0 && ci(2)>0)
            H(ch,tt)=1;
        end
    end
end

Theta_channels=find(sum(H,2)>0);

% 2. Use all of PFC/ channels, specify the channel numbers, or use parcellation values
Pfc_channels=find(ismember(ParcellationValues(:,8),[1:6,8,10:12]));


%% GLM analysis: ERP subtraction

Trials_C=find(~isnan(TrialDet(:,12)) & TrialDet(:,14)==1);
Trials_C0=find(strcmp(StimLocation,'None')==1 & ~isnan(TrialDet(:,12)) & TrialDet(:,14)==1);
Trials_I=find(~isnan(TrialDet(:,12)) & TrialDet(:,14)==2);
Trials_I0=find(strcmp(StimLocation,'None')==1 & ~isnan(TrialDet(:,12)) & TrialDet(:,14)==2);

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

erp=ft_data3.trial{Trials_C0(1)};
for ii=2:length(Trials_C0)
erp=erp+ft_data3.trial{Trials_C0(ii)};
end
ErpC=erp./length(Trials_C0);

for cc=1:length(Trials_C)
   ft_data_C.trial{cc}= ft_data_C.trial{cc}-ErpC;
end

erp=ft_data3.trial{Trials_I0(1)};
for ii=2:length(Trials_I0)
erp=erp+ft_data3.trial{Trials_I0(ii)};
end
ErpI=erp./length(Trials_I0);

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

Id_time0=find(ft_freq.time>-0.6 & ft_freq.time<0);
Id_time1=find(ft_freq.time>0.1 & ft_freq.time<=1.4);

% Use the PFC channels that have enhanced Theta post image onset
channels=intersect(Theta_channels,Pfc_channels);

Theta_C0=nanmean(nanmean(ft_freqC.powspctrm(:,channels,:,Id_time0),4),3);
Theta_C=nanmean(nanmean(ft_freqC.powspctrm(:,channels,:,Id_time1),4),3);

Theta_I0=nanmean(nanmean(ft_freqI.powspctrm(:,channels,:,Id_time0),4),3);
Theta_I=nanmean(nanmean(ft_freqI.powspctrm(:,channels,:,Id_time1),4),3);

save([MainPath,'PSD_',Subject,'_',stimType,'.mat'],'Trials_C','Trials_I','Theta_C','Theta_C0','Theta_I','Theta_I0','channels','Pfc_channels','-append');

%% GLME alalysis: finding optimal stimulation site based on RT.
Trials_NS0=find(strcmp(StimLocation,'None')==1 & ~isnan(TrialDet(:,12)));
BlockStimLoc=unique(BlockStimLocation)';
StimType=[{'R Dorsal'},{'R Ventral'},{'L Dorsal'},{'L Ventral'}];
StimId=[]; 
for ss=1:length(StimType)
    sid=find(strcmp(BlockStimLoc,StimType{ss})==1);
    if(~isempty(sid))
        StimId=[StimId,sid];
    end
end

dRT=zeros(size(StimId));
for ss=1:length(StimId)
     Trials_Stim=find(ismember(StimLocation,BlockStimLoc(StimId(ss)))==1 & ~isnan(TrialDet(:,12)));
     dRT(ss)=nanmean(TrialDet(Trials_NS0,12))-nanmean(TrialDet(Trials_Stim,12));
end
[~,sId]=max(dRT);
Effective_Stim=BlockStimLoc(StimId(sId));

save([MainPath,Subject,'PSD',stimType,'.mat'],'Trials_C','Trials_I','Theta_C','Theta_C0','Theta_I','Theta_I0','channels','Pfc_channels','Effective_Stim','-append');

%% GLME analysis: building the tables

clear
filename=dir('*.mat');
LogTheta=[];
BlockStim=[];
StimAny=[];
Subject=[];
Region=[];
Channel_name=[];

for ff=1:length(filename)
    load(filename(ff).name,'Theta_C','Theta_C0','Theta_I','Theta_I0','StimLocation','MaxTheta_Id','channels','channel_label','Trials_C','Trials_I','TrialDet','Effective_Stim','ParcellationValues');
    
    Trials=[Trials_C;Trials_I];
    NS1=find(strcmp(StimLocation(Trials),'None')==1);
    %NS2=find(strcmp(StimLocation(Trials),'None')==0 & TrialDet(Trials,5)==0);
    %NS2=find(strcmp(StimLocation(Trials),Effective_Stim)==1 & TrialDet(Trials,5)==0);
    NS=[NS1;NS2];
    Stim=[zeros(length(NS1),1);ones(length(NS2),1)];
    
    RThetaC=log(Theta_C./Theta_C0);
    RThetaI=log(Theta_I./Theta_I0);
    RTheta=[RThetaC;RThetaI];
    
    RTheta_NS1=RTheta(NS1,:);
    RTheta_NS2=RTheta(NS2,:);
    RTheta_NS=RTheta(NS,:);
    
    [~,mId2]=max(nanmean(RTheta_NS2));
    [~,mId1]=max(nanmean(RTheta_NS1));
    [~,mId0]=max(nanmean(RTheta_NS));
    
    
    %=find(channels==MaxTheta_Id);
    mId=mId0;
    
    sub=strings(length(NS),1);
    sub(:)=['P',num2str(ff)];
    
    R=strings(length(NS),1);
    R(:)=channel_label(channels(mId));
    Channel_name=[Channel_name;R];
    Region=[Region;ParcellationValues(channels(mId),8)];
    
    LogTheta=[LogTheta;RTheta(NS,mId)];
    BlockStim=[BlockStim;StimLocation(Trials(NS))];
    StimAny=[StimAny;Stim];
    Subject=[Subject;sub];
    
    subplot(3,3,ff)
    plot(nanmean(RTheta_NS1));hold on
    plot(nanmean(RTheta_NS2));box off
    %title([filename(ff).name,' ',Effective_Stim])
    
end

StimType=[{'None'},{'L Ventral'},{'L Dorsal'},{'R Ventral'},{'R Dorsal'}];
StimAny(~ismember(BlockStim,StimType))=-10;

% Fit GLME
% All stim

dataTable = table(LogTheta,categorical(BlockStim,{'None','L Ventral','L Dorsal','R Ventral','R Dorsal'}),categorical(StimAny,[0,1]),Subject,...
         'VariableNames',{'LogTheta','BlockStim','StimAny','Subject'});       

M1 = fitglme(dataTable,'LogTheta ~ BlockStim + (1|Subject)');

% Most effective stim

 dataTable = table(LogTheta,categorical(StimAny),BlockStim,Subject,Region,...
          'VariableNames',{'LogTheta','StimAny','BlockStim','Subject','ChannelName'});       

M2 = fitglme(dataTable,'LogTheta ~ StimAny + (1|Subject)');

%% Theta for first and last blocks
filename=dir('*.mat');
LogTheta =[];
blockNS =[];
sub=[];
    
for ff=1:length(filename)
        load(filename(ff).name,'Theta_C','Theta_C0','Theta_I','Theta_I0','StimLocation','BlockNumber','BlockStimLocation','Trials_C','Trials_I','TrialDet');
        
        Trials=[Trials_C;Trials_I];
        RThetaC=log(Theta_C./Theta_C0);
        RThetaI=log(Theta_I./Theta_I0);
        RTheta=[RThetaC;RThetaI];
        
        NS2=find(strcmp(StimLocation(Trials),'None')==0 & TrialDet(Trials,5)==0);
        RTheta_NS2=RTheta(NS2,:);
        [~,mId]=max(nanmean(RTheta_NS2));
        
        B1=find(strcmp(StimLocation(Trials),'None')==1 & BlockNumber(Trials)==1);
        RTheta_B1=RTheta(B1,mId);
        
        if(strcmp(BlockStimLocation{end},'None')==1)
            B2=find(strcmp(StimLocation(Trials),'None')==1 & BlockNumber(Trials)==max(BlockNumber));
            RTheta_B2=RTheta(B2,mId);
        else
            RTheta_B2=[];
        end
        
        LogTheta=[LogTheta;RTheta_B1;RTheta_B2];
        
        blockNS=[blockNS;zeros(length(RTheta_B1),1);ones(length(RTheta_B2),1)];
        
        S=strings(length(RTheta_B1)+length(RTheta_B2),1);
        S(:)=['P',num2str(ff)];
        sub=[sub;S];
               
end

dataTable = table(LogTheta,categorical(blockNS),sub,...
         'VariableNames',{'LogTheta','BlockNS','Subject'});       

M = fitglme(dataTable,'LogTheta ~ BlockNS + (1|Subject)');

 TH0=LogTheta(blockNS==0);
 TH1=LogTheta(blockNS==1);
 
 figure;
bplot(TH0,1,'std','nomean','width',0.5);hold on
bplot(TH1,2,'std','nomean','width',0.5);
set(gca,'xtick',1:2,'xticklabel',[{'First Block'},{'Last Block'}])

%% Plotting figure 2D
sub=unique(Subject);
P=nan(length(StimType),length(sub));
for ss=1:length(StimType)
    %P=[];
    for pp=1:length(sub)
        P(ss,pp)=nanmean(LogTheta(strcmp(BlockStim,StimType{ss})==1 & strcmp(Subject,sub{pp})==1));
    end
    bplot(P(ss,:),ss,'std','nomean','width',0.6,'box',25);hold on;%scatter(ss.*ones(1,length(P)),P,'k','filled')
end
set(gca,'xtick',1:5,'xticklabel',StimType)
box off
ylabel('Log Theta Ratio')

sub=unique(Subject);
P0=[];P1=[];
for ss=1:length(sub)
    P0(ss)=nanmean(LogTheta(StimAny==0 & strcmp(Subject,sub{ss})==1));
    P1(ss)=nanmean(LogTheta(StimAny==1 & strcmp(Subject,sub{ss})==1));
end
bplot(P0,1,'std','nomean','width',0.6,'box',25);hold on
bplot(P1,2,'std','nomean','width',0.6,'box',25);
set(gca,'xtick',1:2,'xticklabel',[{'None'},{'Capsular Stim'}])
box off
ylabel('Log Theta Ratio')

%% Plotting (figure 2E)
StimLoc=Effective_Stim; % Specify the stimulation to be plotted
Id_freq=find(ft_freq.freq<9 & ft_freq.freq>3);
Id_time0=find(ft_freq.time>-0.6 & ft_freq.time<0);

Trials_NS1=find(strcmp(StimLocation,'None')==1 & ~isnan(TrialDet(:,12)));
Trials_Stim=find(strcmp(StimLocation,StimLoc)==1 & ~isnan(TrialDet(:,12)));
Trials_NS2=intersect(Trials_Stim,find(TrialDet(:,5)==0));
Trials_SS=intersect(Trials_Stim,find(TrialDet(:,5)==1));

%ChId=find(ismember(Pfc_channels,Theta_Stim{ss}));
ChId=MaxTheta_Id;

Theta=squeeze(nanmean(ft_freq.powspctrm(:,ChId,Id_freq,:),3));
Theta_NS1=Theta(Trials_NS1,:);
Theta_NS2=Theta(Trials_NS2,:);
Theta_SS=Theta(Trials_SS,:);

Base_NS1=nanmean(Theta_NS1(:,Id_time0),2);
Base_NS2=nanmean(Theta_NS2(:,Id_time0),2);
Base_SS=nanmean(Theta_SS(:,Id_time0),2);

Time=ft_freq.time;
Id_time1=find(Time>0 & Time<=2);

    P1=[];
    P3=[];P2=[];
    
    for tt=1:length(Id_time1)
        P1(:,tt)=Theta_NS1(:,Id_time1(tt))./Base_NS1;
        P2(:,tt)=Theta_NS2(:,Id_time1(tt))./Base_NS2;
        P3(:,tt)=Theta_SS(:,Id_time1(tt))./Base_SS;
    end
    
    figure(1)
    plot(Time(Id_time1),nanmean(P1),'k','linewidth',2);hold on
    plot(Time(Id_time1),nanmean(P2),'b','linewidth',2);
    plot(Time(Id_time1),nanmean(P3),'r','linewidth',2);box off
    ylabel(channel_label(Theta_Stim{ss}(ch)))

    figure(2)
    plot(Time(Id_time1),nanmean(Theta_NS1(:,Id_time1)),'k','linewidth',2);hold on
    plot(Time(Id_time1),nanmean(Theta_NS2(:,Id_time1)),'b','linewidth',2);
    plot(Time(Id_time1),nanmean(Theta_SS(:,Id_time1)),'r','linewidth',2);box off
    
    %% Theta for diff stim sites
    StimType=[{'R Dorsal'},{'R Ventral'},{'L Dorsal'},{'L Ventral'}];
    clr=['g','m','b','r'];
    
for ff=1:length(filename)
        load(filename(ff).name,'Theta_C','Theta_C0','Theta_I','Theta_I0','StimLocation','channels','channel_label','Trials_C','Trials_I','TrialDet');
        
        Trials=[Trials_C;Trials_I];
        NS1=find(strcmp(StimLocation(Trials),'None')==1);
        RThetaC=log(Theta_C./Theta_C0);
        RThetaI=log(Theta_I./Theta_I0);
        RTheta=[RThetaC;RThetaI];
        RTheta_NS1=RTheta(NS1,:);
        subplot(3,3,ff)
        plot(nanmean(RTheta_NS1),'k','linewidth',1);hold on
        
        for ss=1:length(StimType)
            
            NS2=find(strcmp(StimLocation(Trials),StimType{ss})==1 & TrialDet(Trials,5)==0);
            
            RTheta_NS2=RTheta(NS2,:);
            
            
            plot(nanmean(RTheta_NS2),clr(ss),'linewidth',1);box off
            axis tight
            title(filename(ff).name(1:5))
        end
        
end
legend([{'None'},StimType])

 %% Percentage of channels with increase and decrease in Theta for diff stim sites
 filename=dir('*.mat');
 
    StimType=[{'R Dorsal'},{'R Ventral'},{'L Dorsal'},{'L Ventral'}];
    clr=['g','m','b','r'];
    
   Ntotal=nan(1,length(filename));
   Ninc=nan(length(StimType),length(filename));
   Ndec=Ninc;
   RegionInc=cell(length(StimType),length(filename));
   RegionDec=RegionInc;
   
for ff=1:length(filename)
        load(filename(ff).name,'Theta_C','Theta_C0','Theta_I','Theta_I0','StimLocation','channels','channel_label','Trials_C','Trials_I','TrialDet','ParcellationValues');
        
        Trials=[Trials_C;Trials_I];
        NS1=find(strcmp(StimLocation(Trials),'None')==1);
        RThetaC=log(Theta_C./Theta_C0);
        RThetaI=log(Theta_I./Theta_I0);
        RTheta=[RThetaC;RThetaI];
        RTheta_NS1=RTheta(NS1,:);
        Ntotal(ff)=length(channels);
        
        for ss=1:length(StimType)
            
            NS2=find(strcmp(StimLocation(Trials),StimType{ss})==1 & TrialDet(Trials,5)==0);
            
            RTheta_NS2=RTheta(NS2,:);
            
            Diff_NS2_NS1=nanmean(RTheta_NS2)-nanmean(RTheta_NS1);
            Ninc(ss,ff)=length(find(Diff_NS2_NS1>0));
            Ndec(ss,ff)=length(find(Diff_NS2_NS1<=0));
            
            RegionInc{ss,ff}=ParcellationValues(channels(Diff_NS2_NS1>0),8);
            RegionDec{ss,ff}=ParcellationValues(channels(Diff_NS2_NS1<=0),8);
        end
        
end

Rinc=Ndec./repmat(Ntotal,4,1);
Rdec=Ndec./repmat(Ntotal,4,1);
Rinc(Rinc==0)=nan;
Rdec(Rdec==0)=nan;

RLabel=[{'DLPFC'},{'DLPFCPost'}, {'VLPFC'}, {'DMPFC'}, {'MedOFC'},{'LatOFC'}, {'DACC'}, {'RACC'}, {'PACC'}];
Raxis=[1,3,6,2,4,5,11,10,12];
for ss=1:length(StimType)
    Ri=cell2mat(RegionInc(ss,:)');
    Ri=Ri(~isnan(Ri));Ri(Ri==8)=6;
    Ni=[];Nd=[];
    for rr=1:length(Raxis)
        Ni(rr)=length(find(Ri==Raxis(rr)));
    end
    
    Rd=cell2mat(RegionDec(ss,:)');
    Rd=Rd(~isnan(Rd));Rd(Rd==8)=6;
    for rr=1:length(Raxis)
        Nd(rr)=length(find(Rd==Raxis(rr)));
    end
    
    subplot(2,2,ss)
    bar(Ni,'b');hold on; bar(-Nd,'k')
    set(gca,'xtick',1:9,'xticklabel',RLabel,'xticklabelrotation',25)
    box off
    title(StimType{ss})
end

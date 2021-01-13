clear all

%% This code creates neural features and state estimates for each subject dataset and saves them for encoder-decoder modeling 
MainPath='D:\BasuDDrive\MSIT paper\Neural Data\'; % Specify Data Path
Subject='P01'; % Specify subject code

%% Convert input LFP trials saved in Fieldtrip format  to Power Spectrum. This requires the Fieldtrip toolbox

%stimType='No_Stim2';
stimType='STIM';
filename=['FT data\Trial LFP\BIPOLFieldTripFormat_AlignedToImagePresent_',Subject,'MSIT',stimType,'.mat'];
fileInterictal=['interictal channels\',Subject,'Channels_InterictalActivity_NS2.mat'];

load([MainPath,filename]);
load([MainPath,fileInterictal]);

ChannelPairNames = [ChannelPairNamesBank1 ; ChannelPairNamesBank2];
for c=1:length(ft_data3.label)
   chIdx = str2num(ft_data3.label{c});
   ft_data3.label{c} = ChannelPairNames(chIdx,:);
end

% Get rid of sampleinfo if it exists.
ft_data3 =  rmfield(ft_data3,'sampleinfo');

% Notch filter the trials.
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = [115 135];
cfg.bsfiltord = 4;

addpath 'D:\C Drive\matlab codes\fieldtrip-master'
ft_data3 = ft_preprocessing(cfg,ft_data3);

% PSD calculation
cfg = [];
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.method = 'wavelet';
cfg.foi = [2:50 55:5:200];
cfg.toi = [-1.9:0.01:3.9];
ft_freq = ft_freqanalysis(cfg,ft_data3);

channel_label=ft_freq.label;

savePath=[MainPath,'PSD\'];

save([savePath,Subject,'PSD',stimType,'.mat'],'ft_data3','ft_freq','ch_ictal','channel_label','TrialDet','-v7.3');

% Run if STIM data
BlockNumber = BlockNumber';
BlockType = BlockType';

%% Uncomment the specific P#
%BlockStimLocation = {'None','None','LDacc1','LDacc1','LDacc2','LDacc2','None'}; %P02 CL
%BlockStimLocation = {'None','LVentral','RVentral','LMiddle','RMiddle','None'}; %P08
%BlockStimLocation = {'None','R Ventral','L Dorsal','BAD','R Dorsal','LVentral','None'};%P09
%BlockStimLocation = {'None','L Dorsal','L Ventral','L Dorsal'}; %P10
%BlockStimLocation =  {'None','None','Bad','Bad','Bad','R Ventral','RDorsal','L Ventral','L Middle','L Dorsal'}; %P11
%BlockStimLocation =  {'None','None','R Ventral','R Dorsal','L Ventral','LDorsal','R Dorsal Fix','R Dorsal delay','R Dorsal','None'};% P12
%BlockStimLocation =  {'None','None','R Ventral','R Dorsal','L Ventral','L
%Dorsal','L Dorsal Fix','L Dorsal delay','R Dorsal Fix','R Dorsaldelay','None'}; % P13
%BlockStimLocation ={'None','None','None','R Dorsal','L Dorsal','L Dorsal'};% P14, part1
%BlockStimLocation ={'None','R Ventral','R Dorsal','R Dorsal'};% P14, part2
%BlockStimLocation ={'None','None','aborted','single pulse','R Dorsal','RDorsal','R Ventral','None'}; % P15
%BlockStimLocation={'None','RDorsal','RDorsal','LDorsal','Rventral','RDorsal','None','RDorsal_HGP','RVentral_HGP','None'};%P17

StimLocation = cell(size(TrialDet,1),1);
for b=1:length(BlockStimLocation)
   StimLocation(BlockNumber==b) = deal(BlockStimLocation(b)); 
end
save([savePath,Subject,'PSD',stimType,'.mat'],'ft_data3','ft_freq','ch_ictal','channel_label','TrialDet','BlockStimLocation','StimLocation','BlockNumber','BlockType','-v7.3');

%% Neural feature matrix and behavior

clear all
Subject='P01';
MainPath='D:\BasuDDrive\MSIT paper\Neural Data\';

files=ls(['PSD\',Subject,'*.mat']);

% Behavior
RT=[];
I=[];
seq=[];
seq_id=[];

Nnan=50;
id_rt=12;
id_conflict=14;

for ff=1:size(files,1)
    load([MainPath,'PSD\',files(ff,:)],'ch_ictal','TrialDet');
    RT=[RT;TrialDet(:,id_rt);nan(Nnan,1)];
    seq=[seq,length(TrialDet(:,id_rt)),Nnan];
    seq_id=[seq_id,1,0];
    I=[I;TrialDet(:,id_conflict);zeros(Nnan,1)];
end

RT=RT(1:length(RT)-Nnan);
I=I(1:length(I)-Nnan);
seq=seq(1:end-1);
seq_id=seq_id(1:end-1);

figure;plot(RT)
% Reject channels with interictal activity and find common channels 
channels_exclude=[];
load([MainPath,'PSD\',files(1,:)],'channel_label');
channels_common=channel_label;

for fl=1:size(files,1)
    
    load([MainPath,'PSD\',files(fl,:)],'channel_label','ch_ictal');
    channels_common=intersect(channels_common,channel_label);
    channels_exclude=union(channels_exclude,mat2cell(ch_ictal,ones(1,size(ch_ictal,1)),size(ch_ictal,2)));
end
channel_name=setdiff(channels_common,channels_exclude);

% Neural features finally
F=[4,8;8,15;15,30;30,55;70,110;130,200]; % Frequency bands
T0=0;
Twin=2; % Time window
M=[];   
S=[];

for fl=1:size(files,1)
    
    load([MainPath,'PSD\',files(fl,:)],'ft_freq','TrialDet');
      
    trIdx=1:length(TrialDet);
    
    plotChans=channel_name;
    L=length(plotChans);
    Mean_feature=zeros(length(trIdx),(size(F,1)*length(T0)),L);
    Std_feature=zeros(length(trIdx),(size(F,1)*length(T0)),L);
    
    
    for ch=1:L
        chid=find(strcmp(ft_freq.label,plotChans{ch})==1);
        for ff=1:size(F,1)
            fIdx = find(ft_freq.freq >= F(ff,1) & ft_freq.freq <= F(ff,2));
            
            tIdx = find(ft_freq.time >= T0(1) & ft_freq.time <= T0(1)+Twin);
            powPerTrial =  squeeze(mean(log10(ft_freq.powspctrm(trIdx,chid,fIdx,tIdx)),3));
            Mean_feature(:,ff,ch)= nanmean(powPerTrial');
            Std_feature(:,ff,ch)=  nanstd(powPerTrial');
           
        end
        
    end
  M=[M;Mean_feature];
  S=[S;Std_feature];
end

% Creating feature matrix for model fitting
Y=[];Ys=[];YL=[];

for ff=1:size(F,1)
Y=[Y,squeeze(M(:,ff,:))];
Ys=[Ys,squeeze(S(:,ff,:))];
for ch=1:L
YL=[YL;[ff,channel_name(ch)]];
end
end

%% State estimate from behavior

N = length(RT);
% Yn - log of reaction time
Yn = log(RT);
% Yb - correct/incorrect not used here
Yb = ones(N,1);
% Input, 1 xi
In = zeros(N,2);
In(:,1)=1;
In(find(I==max(I)),2)=1;
% Input, Ib equal to In
Ib = In;
% Uk, which is zero
Uk = zeros(N,1);
% Valid, which is valid for observed point
Valid = zeros(N,1);
Valid(find(isfinite(RT)))=1;

Param = compass_create_state_space(2,1,2,2,eye(2,2),[1 2],[0 0],[1 2],[0 0]);
% set learning parameters
Iter  = 1000;
Param = compass_set_learning_param(Param,Iter,0,1,1,0,1,1,1,2,1);

%% Format the Data
[XSmt,SSmt,Param,XPos,SPos,ML,YP,~]=compass_em([1 0],Uk,In,Ib,Yn,Yb,Param,Valid);

% Plotting to make sure the algorithm converged and everything looks ok
ml=[];
for i=1:Iter
    ml(i)=ML{i}.Total;
end
plot(ml,'LineWidth',2);
ylabel('ML')
xlabel('Iter');

K  = length(XPos);
xm = zeros(K,2);
xb = zeros(K,2);
for i=1:K
    temp=XPos{i};xm(i,1)=temp(1);xm(i,2)=temp(2);
    temp=SPos{i};xb(i,1)=temp(1,1);xb(i,2)=temp(2,2);
end
subplot(2,1,1);ay_plot_bound(1,(1:K),xm(:,1),(xm(:,1)-2*sqrt(xb(:,1)))',(xm(:,1)+2*sqrt(xb(:,1)))');hold on;plot(Yn)
subplot(2,1,2);ay_plot_bound(1,(1:K),xm(:,2),(xm(:,2)-2*sqrt(xb(:,2)))',(xm(:,2)+2*sqrt(xb(:,2)))');

%% Saving data for modeling
stimType='NoStim_Stim';
save([MainPath,'Extracted Features\',Subject,'_features_',stimType,'.mat'],'Y','Ys','YL','seq','seq_id','RT','In','channel_name','Param','XPos','SPos','XSmt','SSmt','YP');

%% Getting parcellations
clear all
load('D:\BasuDDrive\MSIT paper\Neural Data\FT data\parcellation\MG120_NS1.mat')
Reference_channels=[];
ChannelPairs=[ChannelPairNamesBank1;ChannelPairNamesBank2];
for ii=1:length(ChannelPairs)
Reference_channels{ii}=ChannelPairs(ii,:);
end

load('P01_features_No_Stim.mat', 'channel_name')
[~,id1,id2]=intersect(Reference_channels,channel_name);

Regions=ParcellationValues(id1,8);



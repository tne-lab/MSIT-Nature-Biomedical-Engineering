%% This script is to fir an encoder model and then use a feature reduction procedure to determine a reduced model
clear all

ReadFile='D:\BasuDDrive\MSIT paper\Neural Data\Extracted Features\';
SaveFile='D:\BasuDDrive\MSIT paper\Neural Data\Models\';

subject='P09';
stimType='Stim';

file_name = [ReadFile,subject,'_features_',stimType,'.mat'];

%% Training Phase 
x_min = -2;
x_max =  2;
sample= 2000;
Xs    = linspace(x_min,x_max,sample);

ModelSetting.pName             = [SaveFile,subject,'_Model_Trainon_NoStim'];
%%--------------------------------------------------------------
% call this on learning data
% define file name and state variable that being estimated
ModelSetting.pVal             = 0.01;       % 0.05, 0.01, 0.001
ModelSetting.SelMode          = 6;          % 6 or 7
ModelSetting.NoStateSamples   = 1000;       % No of trajectories 
ModelSetting.which_state      = 2;          % 1: Baseline, 2: Conflict
ModelSetting.Xs = Xs;  

% Load file containing neural features and state values
load(file_name);
temp = cell2mat(XPos');

no_feature = size(Y,2);

XM   = temp(ModelSetting.which_state,:);
for i=1:length(XM)
    temp=SPos{i};XS(i)=temp(1,1);
end

% Visualize states to determine training set
plot(XM);

L=length(XM);
Trial_seq=cumsum(seq(find(seq_id==1)));
Lt=round(L/2); % Length of training set, we used L/2 or 2L/3
Start_id=0; % starting index of training set. 
% Other options
%Start_id=Trial_seq(1)0;
%Start_id=Trial_seq(1)+round(L/2)0;

% Assign training and test sets, make sure that the training set is >1
TrainInd=Start_id+1:Start_id+Lt;
TestInd=setdiff(1:L,TrainInd);

% If we want to discard stimulated trials in the training set
NoStimTrials=find(TrialDet(:,5)==0);
TrainInd=intersect(TrainInd,NoStimTrials');

% Specify path where the encoder modeling codes are stored. A compessed
% folder is included in the repository.
% Corresponding paper: https://www.mitpressjournals.org/doi/full/10.1162/neco_a_01196

addpath 'D:\MSIT Backup\GH\Decoder_Encoder_Model_MultipleTrajectory'

ModelName=ay_neural_encoder_training(file_name,ModelSetting,TrainInd);

%%----------------------------------------------------
% load the model
load(ModelName);
ind=find(dValid(:,1)==1);  % If there are too many features that were passed by f-test , can reduce pvalue and run again
XM=XM(ValidId)';
XS=XS(ValidId)';

XMh=XM+2.*sqrt(XS);
XMl=XM-2.*sqrt(XS);

%%---------------------------------------------
% build state-transition distribution
TransP = ones(length(Xs),length(Xs));
for i=1:length(Xs)
    TransP(i,:)=pdf('normal',Xs(i),sParam.a*Xs,sqrt(sParam.sv));
end

%%-----------------------------------------------

figure(1) % Plotting R^2 values of all neural features
plot(dValid(:,4),'LineWidth',2);
hold on
ind =find(dValid(:,1));
plot(ind,dValid(ind,4),'o','LineWidth',2);
box off
title('R^2 and Valid Features')
xlabel('Feature Index')
ylabel('R^2')

figure(2) % Example of neural feature that individually encodes state
[~,m_ind] = max(dValid(:,4));
[~,c_ind] = min(dValid(:,4));
tY = mean(eParam{m_ind}.Y);
fY = Y(TrainInd,m_ind);
subplot(3,1,1)
plot(fY,'LineWidth',2);
hold on
plot(tY,'LineWidth',2);
box off
title(['model ' eParam{m_ind}.RefModel  ', slope(x)=' num2str(eParam{m_ind}.W(2))])
legend('Feature','Prediction')
axis tight
%ylabel(channel_name(Yl(m_ind,2)))
subplot(3,1,2)
plot(TrainInd,XM(TrainInd),'LineWidth',2);
xlabel('Training Index')
ylabel('X');box off
axis tight
subplot(3,1,3)
plot(XM(TrainInd),Y(TrainInd,m_ind),'*');hold on;plot(XM(TrainInd),Y(TrainInd,c_ind),'o');
xlabel('X')
%ylabel(channel_name(Yl(m_ind,2))); 
ylabel('Z')
box off
axis tight
title('Scatter Plot')

%%------------------------------------
%This is the 2nd step for shrinking neural features, using a pruning method

TProb = ay_individual_decoder(data_type,eParam,Xs,dValid(:,1),Y);
XProb = TProb;

% Subset feature given Training Data
for f=1:length(TProb)
    if  TProb{f}.valid
        TProb{f}.prb=TProb{f}.prb(TrainInd,:);
    end
end
[rmse_ind,rmse_curve,optim_curve,winner_list] = ay_sort_decoder_sub(TProb,Xs,dValid(:,1),SampleX(:,TrainInd));
figure(3)
plot(rmse_curve,'LineWidth',2);

[~,ind]=min(rmse_curve);
opt_train=rmse_ind{ind};

%%----------------------------------------------
% Subset feature given Test Data
for f=1:length(XProb)
    if  XProb{f}.valid
        XProb{f}.prb=XProb{f}.prb(TestInd,:);
    end
end
[xrmse_ind,xrmse_curve,xoptim_curve,xwinner_list] = ay_sort_decoder_sub(XProb,Xs,dValid(:,1),SampleX(:,TestInd));
figure(4)
plot(xrmse_curve,'LineWidth',2)

[~,ind]=min(xrmse_curve);
opt_test=xrmse_ind{ind};


% Visualizing what the mean decoded state looks like with the reduced
% features v all features 

% With Training data
opt_id=opt_train; 
tdValid=zeros(length(dValid(:,1)),1);
tdValid(opt_id)=1;
%%---------------------------------------------------
% we keep previous posterior - initialize XPre
XPre = ones(1,size(TransP,1));
%XPre = pdf('normal',Xs,XPos{1}(1),10.*sqrt(SPos{1}(1,1))); 
% we might use previous point of feature
Yprv = zeros(1,no_feature);
% this is the hypothetical real-time loop
% I am keeping mean of estimate here
MEAN=nan(length(XM),1); LOW=nan(length(XM),1); HI=nan(length(XM),1); 
lMap=zeros(length(XM),length(Xs)); fMap=lMap;
for n=1:size(Y,1)
    % load Yk
    Yk = Y(n,:);
    % decoder model 
    % Using subset of features
    if (~isnan(Yk))
    [XPos,CurEstimate,Xll] = ay_one_step_decoder(data_type,eParam,XPre,TransP,Xs,tdValid,Yk,Yprv);
    lMap(n,:) = Xll;
    fMap(n,:) = XPos;
    % next step
    XPre = XPos;
    Yprv = Yk;
    % result
    MEAN(n) = CurEstimate.Mean;
    LOW(n)  = CurEstimate.Bound(1);
    HI(n)   = CurEstimate.Bound(2);
    else
        continue
    end
end
XM_train=MEAN;
XS_train=[LOW;HI];

% With testing data
opt_id=opt_test; 
tdValid=zeros(length(dValid(:,1)),1);
tdValid(opt_id)=1;
%%---------------------------------------------------
% we keep previous posterior - initialize XPre
XPre = ones(1,size(TransP,1));
%XPre = pdf('normal',Xs,XPos{1}(1),10.*sqrt(SPos{1}(1,1))); 
% we might use previous point of feature
Yprv = zeros(1,no_feature);
% this is the hypothetical real-time loop
% I am keeping mean of estimate here
MEAN=nan(length(XM),1); LOW=nan(length(XM),1); HI=nan(length(XM),1); 
lMap=zeros(length(XM),length(Xs)); fMap=lMap;
for n=1:size(Y,1)
    % load Yk
    Yk = Y(n,:);
    % decoder model 
    % Using subset of features
    if (~isnan(Yk))
    [XPos,CurEstimate,Xll] = ay_one_step_decoder(data_type,eParam,XPre,TransP,Xs,tdValid,Yk,Yprv);
    lMap(n,:) = Xll;
    fMap(n,:) = XPos;
    % next step
    XPre = XPos;
    Yprv = Yk;
    % result
    MEAN(n) = CurEstimate.Mean;
    LOW(n)  = CurEstimate.Bound(1);
    HI(n)   = CurEstimate.Bound(2);
    else
        continue
    end
end
XM_test=MEAN;
XS_test=[LOW;HI];

% With all features
% we keep previous posterior - initialize XPre
XPre = ones(1,size(TransP,1));
%XPre = pdf('normal',Xs,Param.X0(1),10.*ssqrt(Param.W0(1,1)));
% we might use previous point of feature
Yprv = zeros(1,no_feature);
% this is the hypothetical real-time loop
% I am keeping mean of estimate here
MEAN=nan(length(XM),1); LOW=nan(length(XM),1); HI=nan(length(XM),1); 
lMap=zeros(length(XM),length(Xs)); fMap=lMap;
for n=1:size(Y,1)
    % load Yk
    Yk = Y(n,:);
    % decoder model 
    % Using all features 
    if (~isnan(Yk))
    [XPos,CurEstimate,Xll] = ay_one_step_decoder(data_type,eParam,XPre,TransP,Xs,dValid(:,1),Yk,Yprv);
    lMap(n,:) = Xll;
    fMap(n,:) = XPos;
    % next step
    XPre = XPos;
    Yprv = Yk;
    % result
    MEAN(n) = CurEstimate.Mean;
    LOW(n)  = CurEstimate.Bound(1);
    HI(n)   = CurEstimate.Bound(2);
    else
        continue
    end
end

XM_all=MEAN;
XS_all=[LOW;HI];
%%--------------------------------------------
% plot figure, State Mean plus Mean+/-Std Estimate


figure(5)

plot(XM,'b','LineWidth',2);hold on; 
plot(XM_all,'k','LineWidth',2);
plot(XM_train,'r','LineWidth',2);
plot(XM_test,'g','LineWidth',2);
plot(XMh,'b--');plot(XMl,'b--')

scatter(TrainInd,XM(TrainInd),'k')
box off
legend('behavior estimate','neural estimate with full model','Reduced Model: Training','Reduced Model: test');


% Metrics for decoder performance
Id_metrics=TestInd(~isnan(XM_all(TestInd))); % Trials at which the metrics are calculated

% Correlation bw mean decoded states
cc1=corrcoef(XM(Id_metrics),XM_all(Id_metrics)); 
cc2=corrcoef(XM(Id_metrics),XM_train(Id_metrics));
cc3=corrcoef(XM(Id_metrics),XM_test(Id_metrics));
CC=[cc1(1,2),cc2(1,2),cc3(1,2)];

% RMSE
rmse1=mean((XM(Id_metrics)-XM_all(Id_metrics)).^2)/(max(XM(Id_metrics))-min(XM(Id_metrics))); 
rmse2=mean((XM(Id_metrics)-XM_train(Id_metrics)).^2)/(max(XM(Id_metrics))-min(XM(Id_metrics)));
rmse3=mean((XM(Id_metrics)-XM_test(Id_metrics)).^2)/(max(XM(Id_metrics))-min(XM(Id_metrics)));
RMSE=[rmse1,rmse2,rmse3];

% Ratio of trials at which the neurally decoded state was withing the confidence
% bounds of the behavior decoded state
Nin1=[length(find((XM_all(Id_metrics)<XMh(Id_metrics))&(XM_all(Id_metrics)>XMl(Id_metrics))))];
Nin2=[length(find((XM_train(Id_metrics)<XMh(Id_metrics))&(XM_train(Id_metrics)>XMl(Id_metrics))))];
Nin3=[length(find((XM_test(Id_metrics)<XMh(Id_metrics))&(XM_test(Id_metrics)>XMl(Id_metrics))))];

Hdr=[Nin1,Nin2,Nin3]./length(XM(Id_metrics));

save(ModelName,'XM','XS','XM_all','XS_all','XM_train','XS_train','XM_test','XS_test','CC','RMSE','Hdr','opt_train','opt_test','TrainInd','TestInd','-append')



%% Shuffling state trajectories
clear all
ReadFile='D:\MSIT Backup\MSIT paper\Neural Data\Extracted Features\';
SaveFile='D:\MSIT Backup\MSIT paper\Neural Data\Shuffled Models\';

files=dir([ReadFile,'*.mat']);
addpath 'D:\MSIT Backup\GH\Decoder_Encoder_Model_MultipleTrajectory'

NPerm=500;
NFeatures=zeros(length(files),NPerm);

for ff=1:length(files)
    file_name=[ReadFile,files(ff).name];
    model_file=dir([SaveFile,files(ff).name(1:5),'*.mat']);
    load(model_file.name,'ModelSetting','training_ind');

    NFeatures(ff,:)=ay_neural_encoder_training_shuffledtrials(file_name,ModelSetting,training_ind,NPerm);
end

%% Finding features that encode states
clear all
subject='P08';
filename=ls(['D:\BasuDDrive\MSIT paper\Neural Data\Extracted Features\',subject,'*.mat']);
load(['D:\BasuDDrive\MSIT paper\Neural Data\Extracted Features\',filename],'YL','channel_name','Regions');

features_train=cell(length(opt_train),3);
features_test=cell(length(opt_test),3);

features_train(:,1:2)=YL(opt_train,:);
features_test(:,1:2)=YL(opt_test,:);

for ch=1:length(opt_train)
    ch_id=find(strcmp(channel_name,features_train(ch,2))==1);
    features_train{ch,3}=Regions(ch_id);
end

for ch=1:length(opt_test)
    ch_id=find(strcmp(channel_name,features_test(ch,2))==1);
    features_test{ch,3}=Regions(ch_id);
end

%% verifying Linear relationship between neural features and states 
x=[XM(training_ind),XM(training_ind).^2];
id=opt_train;
for ii=1:length(id)

y=Y(training_ind,id(ii));
[B1,stats1] = robustfit(x(:,1),y);
[B2,stats2] = robustfit(x,y);

subplot(4,4,ii)
scatter(XM(training_ind),Y(training_ind,id(ii)));hold on
z1=B1(1)+B1(2)*x(:,1);
z2=B2(1)+B2(2)*x(:,1)+B2(3)_*x(:,2);
rm1=sqrt(nanmean((y-z1).^2));
rm2=sqrt(nanmean((y-z2).^2));

plot(x(:,1),z1,'k', x(:,1),z2,'m','linewidth',2)
title(['p:',num2str(round(stats1.p(2),3)),', ',num2str(round(stats2.p(3),3))])
%title(['RMSE:',num2str(rm1),', ',num2str(rm2)])
axis tight
end
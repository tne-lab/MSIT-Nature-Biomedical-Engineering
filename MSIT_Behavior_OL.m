%% This script is for analyzing and plotting the effect of open loop capsular stimulation on Behavior (Reaction Times, Baseline and Conflict States)
clear all

% Load Behavior Data
load('Behavior_Stim_Subjects.mat')

% Subjects who performed Open Loop stimulation
allSubj = [{'P8'};{'P9'};{'P10'};{'P11'};{'P12'};{'P13'}];
subject_id=find(ismember(subject,allSubj)==1);

% Use data only for open loop experiment subjects
blockNum=blockNum(subject_id);
blockStim=blockStim(subject_id);
interference=interference(subject_id);
responseCorrect=responseCorrect(subject_id);
responseTimes=responseTimes(subject_id);
StimAny=StimAny(subject_id);
subject=subject(subject_id);
switchType=switchType(subject_id);
trialNum=trialNum(subject_id);
trialStim=trialStim(subject_id);

%% Make data into table for modeling the open loop stimulation dataset

dataTable = table(responseTimes,categorical(interference),categorical(StimAny),...
    blockNum,categorical(blockStim,{'None','L VS','L DS','R VS','R DS'}),...
    categorical(trialStim,{'None','L VS','L DS','R VS','R DS'}),...
    trialNum,categorical(switchType,{'CC','CI','IC','II'}),(1:length(blockNum))',subject, ...
    'VariableNames',{'responseTimes','interference','anyStim','blockNum','blockStim',...
    'trialStim','trialNum','switchType','totalTrialNum','subject'});

%% Fit GLME.
% Check if we should use block level or trial level stimulation coding
M1 = fitglme(dataTable,'responseTimes ~ interference + blockStim + blockNum + (1|subject)','Link','log');
M2 = fitglme(dataTable,'responseTimes ~ interference + trialStim + blockNum + (1|subject)','Link','log');
disp(sprintf('AIC %.3f -> %.3f',M1.ModelCriterion.AIC,M2.ModelCriterion.AIC)); % Check if we should use block level or trial level stimulation coding
M=M1;

% Calculate pvalues
pval=M.Coefficients.pValue;
pval=pval(4:end);
[~,~,pfdr] = fdr(pval);

%% Box Plots for RT (Figure 2C)
responseCopy = responseTimes; % RT

% For All types of Stim Effect
plotStrs = {'None','L Ventral','L Dorsal','R Ventral','R Dorsal'};
findStrs = {'None','L VS','L DS','R VS','R DS'};

plotIdx = {};

SubjectMean=nan(length(allSubj),length(plotStrs));

% Calculating on a Subject level 
for s=1:length(allSubj)
    for s1=1:length(plotStrs)
                id = find(strcmp(subject,allSubj{s})& strcmp(blockStim,findStrs{s1}));
        if(~isempty(id))
        SubjectMean(s,s1) = nanmean(responseCopy(id));
        end
    end
end

for s=1:length(plotStrs)
    plotIdx{s} = find(strcmp(blockStim,findStrs{s}));
    barHeights{s} = responseCopy( plotIdx{s});
end

L=cellfun(@length, barHeights);
RMtx=nan(max(L),length(L));
for ll=1:length(L)
    RMtx(1:length(barHeights{ll}),ll)=barHeights{ll};
end

figure;
bplot(RMtx,'std','nomean','width',0.6)
hold on;

axStrs = plotStrs;
for l=1:length(axStrs)
    axStrs{l} = [axStrs{l} ', n=' num2str(length(unique(subject(cellfun(@(x)~isempty(strfind(x,findStrs{l})),blockStim)))))];
end
set(gca,'xtick',1:5,'xticklabel',axStrs,'xticklabelrotation',15)

for ii=1:length(allSubj)
    scatter(1:5,SubjectMean(ii,:),'k','filled');
end

%% Fit state-space model using COMPASS toolbox.
nIter = 1000;
% Model based on Reaction time only
XPos_collect = zeros(2,length(responseTimes)); 
SPos_collect = zeros(2,length(responseTimes));

% Model based on reaction time and accuracy
XPos_collect2 = zeros(2,length(responseTimes));
SPos_collect2 = zeros(2,length(responseTimes));

%Add path to COMPASS toolbox folder
addpath('C:\Users\basuia\OneDrive - University of Cincinnati\Research\MSIT Nature Biomed Engg\Codes\COMPASS_StateSpaceToolbox')
for s=1:length(allSubj)
    
    % Indexes for this one subject.
    sIdx = find(strcmp(subject,allSubj{s}));
    
    Param = compass_create_state_space(2,1,2,2,eye(2,2),[1 2],[0 0],[1 2],[0 0]);
    Param.Cut_Time = log(2);
    
    % We do try to re-learn parameters from data, but not the behavior model.
    Param = compass_set_learning_param(Param,nIter,0,1,1,0,1,1,1,2,1);
    
    RT=responseTimes(sIdx);
    N = length(RT);
    % Yn - log of reaction time
    Yn = log(RT);
    % Yb - correct/incorrect not used here
    Yb = ones(N,1);
    % Input, 1 xi
    In = zeros(N,2);
    In(:,1)=1;
    In(:,2)=interference(sIdx);
    % Input, Ib equal to In
    Ib = In;
    % Uk, which is zero
    Uk = zeros(N,1);
    
    Valid = zeros(N,1);
    Valid(find(isfinite(RT)))=1;
    
    [XSmt,SSmt,Param,XPos,SPos,ML,YP,~]=compass_em([1 0],Uk,In,Ib,Yn,Yb,Param,Valid);
    % Extract model output.
    for i=1:2
        XPos_collect(i,sIdx) =  cellfun(@(x) x(i),XPos);
        SPos_collect(i,sIdx) =  cellfun(@(x) x(i,i),SPos);
    end
    
    % Including Binary error (Accuracy)
    Yb=responseCorrect(sIdx);
    Param = compass_create_state_space(2,1,2,2,eye(2,2),[1 2],[1 1],[1 2],[0 0]); 
    Param = compass_set_learning_param(Param,nIter,0,1,1,1,1,1,1,2,1);
    [XSmt,SSmt,Param,XPos,SPos,ML2,YP,~]=compass_em([1 1],Uk,In,Ib,Yn,Yb,Param,Valid);
    for i=1:2
        XPos_collect2(i,sIdx) =  cellfun(@(x) x(i),XPos);
        SPos_collect2(i,sIdx) =  cellfun(@(x) x(i,i),SPos);
    end
    
end

%% Box Plots for the state variables. Here we use a permutation test for calculating p-values (Figures 3 C,D)

responseCopy = XPos_collect(1,:); % Baseline state
%responseCopy = XPos_collect(2,:); % Conflict State

% For All types of Stim Effect
plotStrs = {'None','L Ventral','L Dorsal','R Ventral','R Dorsal'};
findStrs = {'None','L VS','L DS','R VS','R DS'};

plotIdx = {};

SubjectMean=nan(length(allSubj),length(plotStrs));

% Normalizing within subject
for s=1:length(allSubj)
    sIdx = find(strcmp(subject,allSubj{s}));
    responseCopy(sIdx) = (responseCopy(sIdx) - mean(responseCopy(sIdx)))./std(responseCopy(sIdx));
   
end

% Calculating on a Subject level 
for s=1:length(allSubj)
    for s1=1:length(plotStrs)
                id = find(strcmp(subject,allSubj{s})& strcmp(blockStim,findStrs{s1}));
        if(~isempty(id))
        SubjectMean(s,s1) = nanmean(responseCopy(id));
        end
    end
end


for s=1:length(plotStrs)
    
    plotIdx{s} = find(strcmp(blockStim,findStrs{s}));
    barHeights{s} = responseCopy( plotIdx{s});

    for s2=1:length(plotStrs)
       compIdx = find(strcmp(blockStim,findStrs{s2}));
       [barProbs(s,s2),~,~] = permutationTest(responseCopy(plotIdx{s}),responseCopy(compIdx),1000);
    end

end

% Calculate p-value and fdr correct them
[~,~,pfdr] = fdr(squeeze(barProbs(1,2:end)));

L=cellfun(@length, barHeights);
RMtx=nan(max(L),length(L));
for ll=1:length(L)
    RMtx(1:length(barHeights{ll}),ll)=barHeights{ll};
end

figure;
bplot(RMtx,'std','nomean','width',0.6)
hold on;

axStrs = plotStrs;
for l=1:length(axStrs)
    axStrs{l} = [axStrs{l} ', n=' num2str(length(unique(subject(cellfun(@(x)~isempty(strfind(x,findStrs{l})),blockStim)))))];
end
set(gca,'xtick',1:5,'xticklabel',axStrs,'xticklabelrotation',15)

for ii=1:length(allSubj)
    scatter(1:5,SubjectMean(ii,:),'k','filled');
end

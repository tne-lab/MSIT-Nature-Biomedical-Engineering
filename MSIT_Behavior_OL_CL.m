%% This script is for analyzing and comparing the effect of open loop v closed loop capsular stimulation on Behavior (Reaction Times, Baseline and Conflict States)
clear all

load('Behavior_Stim_Subjects.mat');

%% Fit state-space model.

nIter = 1000;
XPos_collect = zeros(2,length(responseTimes));
SPos_collect = zeros(2,length(responseTimes));

allSubj = unique(smallBlockSubject);
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
    likelihood{s}=ML;
    % Extract model output.
    for i=1:2
        XPos_collect(i,sIdx) =  cellfun(@(x) x(i),XPos);
        SPos_collect(i,sIdx) =  cellfun(@(x) x(i,i),SPos);
    end
    
end

%% Comparing behavior and states during open loop (OL) v closed loop (CL) MSIT stimulation

allSubj = unique(subject);
findStrs={'L DS','L DS CL','R VS','R VS CL','R DS','R DS CL'};
NormStrs = {'None','NoneCL','None','NoneCL','None','NoneCL'};
SubjectMean=nan(length(allSubj),length(findStrs));
barHeights={};

responseCopy = responseTimes; 

% For I-C, use this for getting a conflict variable from RT
conflict=zeros(size(responseCopy));
for s=1:length(allSubj)
    sIdx = find(strcmp(subject,allSubj{s}));
    BN=unique(blockNum(sIdx));
    for bb=1:length(BN)
        bId=find(blockNum==BN(bb)& strcmp(subject,allSubj{s})==1);
        rt=responseCopy(bId);I=interference(bId);
        rtC = mean(rt(I==0));
        conflict(bId) = responseCopy(bId) - rtC; 
    end
end
responseCopy=conflict;

%responseCopy = XPos_collect(1,:);
%responseCopy = XPos_collect(2,:);

for s=1:length(allSubj)
    sIdx = find(strcmp(subject,allSubj{s}));
    responseCopy(sIdx) = (responseCopy(sIdx) - mean(responseCopy(sIdx)))./std(responseCopy(sIdx));
   
end

% P values from permutation test between OL and CL behavior
barProbs=[];
for s=1:2:length(findStrs)
    plotIdx = find(strcmp(blockStim,findStrs{s})==1);
    Idn= find(strcmp(blockStim,NormStrs{s})==1);
    response=responseCopy( plotIdx);
    
    compIdx = strcmp(blockStim,findStrs{s+1});
    comp_response=responseCopy(compIdx)+1-nanmean(responseCopy(strcmp(blockStim,NormStrs{s+1})));
    [barProbs((s+1)/2),~,~] = permutationTest(response+1-nanmean(responseCopy(Idn)),comp_response,1000);
    
end

[~,~,pfdr]=fdr(barProbs);

% plotting (Figures S3, 4 B,C)
barHeights=[];
for s=1:length(findStrs)
    plotIdx = find(strcmp(blockStim,findStrs{s})==1);
    Idn= find(strcmp(blockStim,NormStrs{s})==1);
    response=responseCopy( plotIdx);
    barHeights{s} = response + (1-nanmean(responseCopy(Idn)));
    for p=1:length(allSubj)
        Ids=find(strcmp(subject(plotIdx),allSubj{p})==1);
        if(~isempty(Ids))
            %SubjectMean(p,s) = nanmean(responseCopy(id1))+ (1-nanmean(responseCopy(Idn)));
            SubjectMean(p,s) = nanmean(barHeights{s}(Ids));
        end
    end
end

xt=[1,1.7,3,3.7,5,5.7];

for ii=1:length(barHeights)
bplot(barHeights{ii},xt(ii),'std','width',0.6)
hold on;
end

for ii=1:length(allSubj)
    scatter(xt,SubjectMean(ii,:),'k','filled');
end

xl={'L Dorsal, n=5,2','R Ventral, n=4,3','R Dorsal, n=4,3'}
set(gca,'xtick',[1.25,3.25,5.25])
set(gca,'xticklabel',xl)
set(gca,'xticklabelrotation',0)

%% This script is for analyzing and comparing the effect of open loop v closed loop capsular stimulation on Behavior (Reaction Times, Baseline and Conflict States)
clear all

%% Import csv files containing behavior. Mg 95-105 performed open loop stimulation, MG 106, 108, 112 performed closed loop to behavior stimulation
dataFiles = {
    'MG95/mg95_msit_no-stim_06_2016_01_25_16_08_45.csv',...% no stim 
    'MG95/mg95_msit_num-stim_01_2016_01_25_16_17_34.csv',... % LVF 2-3
    'MG95/mg95_msit_num-stim_02_2016_01_25_16_27_15.csv',... % RVF 3-4
    'MG95/mg95_msit_num-stim_03_2016_01_25_16_34_57.csv',... % LVF 3-4
    'MG95/mg95_msit_num-stim_04_2016_01_25_16_44_16.csv',... % RVF 4-5
    'MG95/mg95_msit_no-stim_07_2016_01_25_16_52_22.csv',... % No stim
     ...
    'MG96/mg96_msit_no-stim_01_2016_02_25_12_00_24.csv',... % no stim 
    'MG96/mg96_msit_num-stim_01_2016_02_25_12_08_16.csv',... %R VS - RVF2-3
    'MG96/mg96_msit_num-stim_02_2016_02_25_13_02_38.csv',... % L DS - LVF 5-6
    'MG96/mg96_msit_num-stim_03_2016_02_25_13_14_41.csv',... % R DS - RVF5-6
    'MG96/mg96_msit_num-stim_04_2016_02_25_13_22_34.csv',... %L VS - LVF 3-4
    'MG96/mg96_msit_no-stim_02_2016_02_25_13_29_04.csv',... % no stim 
    ...
    'MG99/mg99_msit_no-stim_01_2016_04_11_14_16_12.csv',...  %calibration
    'MG99/mg99_msit_num-stim_01_2016_04_11_14_24_20.csv',... % L dorsal LVF 5-6
    'MG99/mg99_msit_num-stim_02_2016_04_11_14_31_05.csv',... % L ventral LVF 3-4
    'MG99/mg99_msit_num-stim_03_2016_04_11_14_40_51.csv',... % L dorsal, notes say "81-82" bad, but not sure what that is... maybe 17-18.
    'MG99/mg99_msit_no-stim_01_2016_04_11_15_58_34.csv',...  % another calibration
    ...    
    'MG102/mg102_msit_no-stim_64_01_2016_05_24_13_40_26.csv',...    % calibration
    'MG102/mg102_msit_no-stim_32_01_2016_05_24_13_46_28.csv',...    % extended calibration
    'MG102/mg102_msit_num-stim_32_03_2016_05_24_14_14_41.csv',...   % RVF4-5, R ventral but in weird place
    'MG102/mg102_msit_num-stim_32_04_2016_05_24_14_19_47.csv',...   % RVF7-8, R dorsal but in weird place
    'MG102/mg102_msit_num-stim_32_05_2016_05_24_14_26_19.csv',...   % LVF5-6, L ventral
    'MG102/mg102_msit_num-stim_32_07_2016_05_24_14_36_05.csv',...   % LVF8-8, L dorsal
    'MG102/mg102_msit_decider-stim_64_01_2016_05_24_14_55_04.csv',...   % trigger on L cingulate HGP being low
    'MG102/mg102_msit_decider-stim_64_02_2016_05_24_15_26_49.csv',...   % trigger on L medial temporal beta being low, adaptive
    'MG102/mg102_msit_no-stim_32_03_2016_05_24_15_33_38.csv',...        % close-out/recovery block
    ...
    'MG104/mg104_msit_no-stim_64_01_2016_08_04_10_58_33.csv',...    % calibration
    'MG104/mg104_msit_no-stim_32_01_2016_08_04_11_04_29.csv',...    % extended calibration
    'MG104/mg104_msit_num-stim_32_01_2016_08_04_11_11_50.csv',...   % R ventral,  RVF4-5
    'MG104/mg104_msit_num-stim_32_02_2016_08_04_11_16_00.csv',...   % R dorsal, RVF8-9
    'MG104/mg104_msit_num-stim_32_03_2016_08_04_11_20_17.csv',...   % L ventral, LVF2-3
    'MG104/mg104_msit_num-stim_32_04_2016_08_04_11_25_10.csv',...   % L dorsal, LVF3-4 (based on cingulate gamma)
    'MG104/mg104_msit_fix-stim_32_01_2016_08_04_11_29_56.csv',...   % R dorsal fixation-locked
    'MG104/mg104_msit_delay-stim_32_01_2016_08_04_11_33_20.csv',... % R dorsal stim-locked + 400ms
    'MG104/mg104_msit_num-stim_32_05_2016_08_04_11_36_36.csv',...   % R dorsal back to num-locked??? or is this actually no-stim?
    'MG104/mg104_msit_no-stim_64_02_2016_08_04_11_40_09.csv',...
    ...
    'MG105/mg105_msit_no-stim_64_01_2016_08_10_14_13_36.csv',...   No stim
    'MG105/mg105_msit_no-stim_64_02_2016_08_10_14_19_34.csv',...    No stim
    'MG105/mg105_msit_num-stim_32_01_2016_08_10_14_27_53.csv',...   R ventral, RVF3-4
    'MG105/mg105_msit_num-stim_32_02_2016_08_10_14_31_37.csv',...   R dorsal, RVF7-8
    'MG105/mg105_msit_num-stim_32_03_2016_08_10_14_35_20.csv',...   L ventra, LVF2-3
    'MG105/mg105_msit_num-stim_32_04_2016_08_10_14_40_26.csv',...   L dorsal, LVF8-9
    'MG105/mg105_msit_fix-stim_32_01_2016_08_10_14_43_44.csv',...   L dorsal, fixation-triggered
    'MG105/mg105_msit_delay-stim_32_01_2016_08_10_14_46_57.csv',... L dorsal, delayed 400ms
    'MG105/mg105_msit_fix-stim_32_02_2016_08_10_14_50_46.csv',...   R dorsal, fixation-triggered
    'MG105/mg105_msit_delay-stim_32_02_2016_08_10_14_53_43.csv',... R dorsal, delayed 400 ms
    'MG105/mg105_msit_no-stim_64_03_2016_08_10_14_59_02.csv',...    No stim
    ...
    'MG106/mg106_msit_no-stim_64_01_2016_11_16_12_07_43.csv',... % no stim
    'MG106/mg106_msit_no-stim_32_01_2016_11_16_12_13_04.csv',... % no stim
    'MG106/mg106_msit_no-stim_32_02_2016_11_16_12_16_27.csv',... % no stimulation, continued practice block
    'MG106/mg106_msit_no-stim_32_03_2016_11_16_12_21_40.csv',... % RVF6-7, 130 Hz, 2 mA, 600 msec, Draper GUI controlled, threshold -0.6, 
    'MG106/mg106_msit_no-stim_32_04_2016_11_16_12_26_12.csv',... % LVF5-6, 130 Hz, 2 mA, 600 msec, Draper GUI controlled, threshold -0.6
    'MG106/mg106_msit_no-stim_64_02_2016_11_16_12_30_14.csv',...  % LVF5-6, 130 Hz, 2 mA, 600 msec, Draper GUI controlled, threshold -0.7
    'MG106/mg106_msit_no-stim_64_03_2016_11_16_13_52_31.csv',... % no stim
    'MG106/mg106_msit_no-stim_32_01_2016_11_16_13_59_47.csv',... % RVF1-2, 130 Hz, 2 mA, 600 msec, Draper GUI controlled, threshold -0.6, 
    'MG106/mg106_msit_no-stim_32_02_2016_11_16_14_17_37.csv',...  % RVF6-7, 130 Hz, 2 mA, 600 msec, Draper GUI controlled, threshold -0.7,  
    ...
    'MG108/mg108_msit_no-stim_64_01_2017_04_10_17_02_55.csv',...% B1: no stim block1
    'MG108/mg108_msit_no-stim_32_01_2017_04_10_17_08_34.csv',... % B2: no stim
    'MG108/mg108_msit_no-stim_32_03_2017_04_10_17_18_11.csv',... % B5: RVF 7-8, initially not CL
    'MG108/mg108_msit_no-stim_32_04_2017_04_10_17_23_04.csv',... % B6: RVF 7-8
    'MG108/mg108_msit_no-stim_32_05_2017_04_10_17_28_07.csv',...  % B7: RVF 4-5
    'MG108/mg108_msit_no-stim_32_06_2017_04_10_17_36_29.csv',...  % B8: No Stim
    'MG108/mg108_msit_no-stim_32_02_2017_04_10_17_48_40.csv',...  % B10:  RVF 5-6  Decider
   ...
    'MG112/stim/mg106_msit_no-stim_64_01_2017_09_06_11_46_15_01.csv',...% B1: no stim block1
    'MG112/stim/mg106_msit_no-stim_32_01_2017_09_06_11_54_31_02.csv',... % B2: RVF 7-8, threshold -0.08
    'MG112/stim/mg106_msit_no-stim_32_02_2017_09_06_12_04_51_03.csv',... % B3: RVF 7-8, theshold  -0.06 -> -0.04
    'MG112/stim/mg106_msit_no-stim_32_03_2017_09_06_12_14_06_04.csv',... % B4: LVF 7-8, theshold  -0.04
    'MG112/stim/mg106_msit_no-stim_32_04_2017_09_06_12_21_30_05.csv',... % B5: RVF 3-4, theshold  -0.06
    'MG112/stim/mg106_msit_no-stim_32_05_2017_09_06_12_27_02_06.csv',... % B6: RVF 7-8, theshold  -0.06
    'MG112/stim/mg106_msit_no-stim_32_06_2017_09_06_13_13_50_07.csv',...  % B7: No Stim
    'MG112/stim/mg106_msit_no-stim_32_01_2017_09_06_13_19_02_08.csv',...  % B8: RVF 7-8 with decider to dacc HGP
    'MG112/stim/mg106_msit_no-stim_32_03_2017_09_06_13_32_13_09.csv',... % B9: RVF 3-4, decider 
    'MG112/stim/mg106_msit_no-stim_32_02_2017_09_06_13_28_31_10.csv',...  % B10:  No Stim
}


% Block Numbers.
smallBlockNum = [1 2 3 4 5 6,...             % mg95
                 1 2 3 4 5 6,...             % mg96
                 1 2 3 4 5,...               % mg99   % 2 3
                 1 2 3 4 5 6 7 8 9,...       % mg102
                 1 2 3 4 5 6 7 8 9 10,...    % mg104
                 1 2 3 4 5 6 7 8 9 10 11,...  % mg105
                  1 2 3 4 5 6 7 8 9,...        % mg106
                  1 2 3 4 5 6 7 ,...          % mg108
                  1 2 3 4 5 6 7 8 9 10,...       % mg112
                 ]        
                                        
% Block Stim. the closed loop ones are indicated by CL
smallBlockStim = {'None','L VS','R other','L Mid','R other','None',...
                  'None','R VS','L DS','R DS','L VS','None',...
                  'None','L DS','L VS','L DS','None',...  
                  'None','None','R VS','R DS','L VS','L DS','L DS CL Cing','L DS CL Temp','None',...
                  'None','None','R VS','R DS','L VS','L DS','R DS fix','R DS delay','R DS','None',...
                  'None','None','R VS','R DS','L VS','L DS','L DS fix','L DS delay','R DS fix','R DS delay','None',...
                  'NoneCL','NoneCL','NoneCL','R DS CL','L DS CL','L DS CL','NoneCL','R VS CL','R DS CL',...
                  'NoneCL','NoneCL','R DS CL','R DS CL','R VS CL','NoneCL','CLPhys',...
                  'NoneCL','R DS CL','R DS CL','L DS CL','R VS CL','R DS CL','NoneCL','R DS CL','R VS CL','NoneCL',...
                  };   
% Subject names              
smallBlockSubject = {'MG95','MG95','MG95','MG95','MG95','MG95',...
                    'MG96','MG96','MG96','MG96','MG96','MG96',...
                    'MG99','MG99','MG99','MG99','MG99',...
                    'MG102','MG102','MG102','MG102','MG102','MG102','MG102','MG102','MG102',...
                    'MG104','MG104','MG104','MG104','MG104','MG104','MG104','MG104','MG104','MG104',...
                    'MG105','MG105','MG105','MG105','MG105','MG105','MG105','MG105','MG105','MG105','MG105',...
                     'MG106','MG106','MG106','MG106','MG106','MG106','MG106','MG106','MG106',...
                     'MG108','MG108','MG108','MG108','MG108','MG108','MG108',...
                     'MG112','MG112','MG112','MG112','MG112','MG112','MG112','MG112','MG112','MG112',...
                    };
                    

%% Extracting the behavior data.

responseTimes = [];
StimAny = [];
interference = [];
blockNum = [];
blockStim = {};
trialStim = {};
trialNum = [];
responseCorrect = [];
switchType = {};
subject = {};

for f=1:length(dataFiles)

    fData = csvread(dataFiles{f},1);
    
    fHandle = fopen(dataFiles{f});
    fHeader = textscan(fHandle,'%[^0123456789,]','Delimiter',',');
    fHeader = fHeader{1};
    fclose(fHandle);
    
     nTrials = size(fData,1);
     lastIdx = length(responseTimes);
     putIdx = lastIdx + (1:nTrials);
     
     responseTimes(putIdx,1) = fData(:,find(strcmp(fHeader,'ResponseTime')));
     interference(putIdx,1) = fData(:,find(strcmp(fHeader,'Condition'))) == 2;
     StimAny(putIdx,1) = ~strcmp(smallBlockStim{f},'None');
     blockNum(putIdx,1) = smallBlockNum(f);
     blockStim(putIdx,1) = deal(smallBlockStim(f));
     
     for i=1:nTrials
         if (fData(i,find(strcmp(fHeader,'FixStimLength'))) > 0) | ...f
              (fData(i,find(strcmp(fHeader,'NumStimLength'))) > 0) | ...
              (fData(i,find(strcmp(fHeader,'RespStimLength'))) > 0)
             trialStim{putIdx(i),1} = smallBlockStim{f};
         else
             trialStim{putIdx(i),1} = 'None';
         end
     end
     
     switchType(putIdx(fData(:,find(strcmp(fHeader,'Conflict'))) == 0),1) = deal({'CC'});
     switchType(putIdx(fData(:,find(strcmp(fHeader,'Conflict'))) == 1),1) = deal({'IC'});
     switchType(putIdx(fData(:,find(strcmp(fHeader,'Conflict'))) == 2),1) = deal({'CI'});
     switchType(putIdx(fData(:,find(strcmp(fHeader,'Conflict'))) == 3),1) = deal({'II'});
     
     trialNum(putIdx,1) = 1:nTrials;
     responseCorrect(putIdx,1) = fData(:,find(strcmp(fHeader,'ResponseAccuracy')));
     subject(putIdx,1) = deal(smallBlockSubject(f));
end
 
%% Data cleanup.

% Remove incorrect and missed trials
keepTrials = responseTimes > 0 & responseCorrect>0;

% Remove known-bad trials if any.
% MG99 not-paying-attention.
badTrials = strcmp(subject,'MG99') & strcmp(blockStim,'L DS') & blockNum == 3 & ismember(trialNum,[17 18]);
keepTrials = keepTrials & ~badTrials;

keepVars = {'responseTimes','StimAny','interference','blockNum',...
    'blockStim','trialStim','trialNum','responseCorrect','switchType','subject'};
for v=1:length(keepVars)
   eval([keepVars{v} ' = ' keepVars{v} '(keepTrials);']); 
end

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

allSubj = unique(smallBlockSubject);
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

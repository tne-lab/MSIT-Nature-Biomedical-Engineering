%% This script is to evaluate Encoder model fit such as Error Metrics

Subject_Stim=[{'P08'},{'P09'},{'P10'},{'P11'},{'P12'},{'P13'},{'P14'},{'P15'},{'P17'}];
files_base=dir('E:\MSIT data for NatureBME\PSD\Models_Base\*.mat');
files_conf=dir('E:\MSIT data for NatureBME\PSD\Models_Conflict\*.mat');

% Finding number of subjects that have a feature in specific region and freq band

Stim_base=zeros(17,6);
NStim_base=zeros(17,6);

Stim_conf=zeros(17,6);
NStim_conf=zeros(17,6);

for ff=1:length(files_base)
    load(files_base(ff).name,'features_test');
    freq=unique(cell2mat(features_test(:,1)));
    reg=unique(cell2mat(features_test(:,3)));
    reg=reg(~isnan(reg));
    if(contains(files_base(ff).name,Subject_Stim))
        Stim_base(reg,freq)=Stim_base(reg,freq)+1;
    else
        NStim_base(reg,freq)=NStim_base(reg,freq)+1;
    end
end
Stim_base(6,:)=Stim_base(6,:)+Stim_base(8,:);
Stim_base=Stim_base([1:7,9:17],:);

NStim_base(6,:)=NStim_base(6,:)+NStim_base(8,:);
NStim_base=NStim_base([1:7,9:17],:);

subplot(1,2,1);imagesc(NStim_base);set(gca,'ytick',1:16,'yticklabel',RegionLabel,'xtick',1:6,'xticklabel',FreqLabel,'xticklabelrotation',25);caxis([0 7])
subplot(1,2,2);imagesc(Stim_base);set(gca,'ytick',[],'xtick',1:6,'xticklabel',FreqLabel,'xticklabelrotation',25);caxis([0 7]);colorbar

for ff=1:length(files_conf)
    load(files_conf(ff).name,'features_test');
    freq=unique(cell2mat(features_test(:,1)));
    reg=unique(cell2mat(features_test(:,3)));
    reg=reg(~isnan(reg));
    if(contains(files_base(ff).name,Subject_Stim))
        Stim_conf(reg,freq)=Stim_conf(reg,freq)+1;
    else
        NStim_conf(reg,freq)=NStim_conf(reg,freq)+1;
    end
end
Stim_conf(6,:)=Stim_conf(6,:)+Stim_conf(8,:);
Stim_conf=Stim_conf([1:7,9:17],:);

NStim_conf(6,:)=NStim_conf(6,:)+NStim_conf(8,:);
NStim_conf=NStim_conf([1:7,9:17],:);

subplot(1,2,1);imagesc(NStim_conf);set(gca,'ytick',1:16,'yticklabel',RegionLabel,'xtick',1:6,'xticklabel',FreqLabel,'xticklabelrotation',25);caxis([0 7])
subplot(1,2,2);imagesc(Stim_conf);set(gca,'ytick',[],'xtick',1:6,'xticklabel',FreqLabel,'xticklabelrotation',25);caxis([0 7]);colorbar

% Collecting performance metrics, features

Header=[{'Stim status'},{'HDR'},{'Correlation coeff'},{'RMSE'},{'No of features'}];
Metric_base=zeros(length(files_base),5);
Metric_conflict=zeros(length(files_base),5);

for ff=1:length(files_base)
    load(files_base(ff).name,'opt_test','CC','Hdr','RMSE');
    
    if(contains(files_base(ff).name,Subject_Stim))
        Metric_base(ff,1)=1;
    end
    
    Metric_base(ff,2:5)=[Hdr,CC,RMSE,length(opt_test)];
end

for ff=1:length(files_conf)
    load(files_conf(ff).name,'opt_test','CC','Hdr','RMSE');
    
    if(contains(files_conf(ff).name,Subject_Stim))
        Metric_conflict(ff,1)=1;
    end
    
    Metric_conflict(ff,2:5)=[Hdr,CC,RMSE,length(opt_test)];
end

% Boxplots of No of Features
ids=find(Metric_base(:,1)==1);idns=find(Metric_base(:,1)==0);
bplot(Metric_base(:,5),1,'width',0.5);hold on; 
scatter(ones(length(ids),1),Metric_base(ids,5),'r*');scatter(ones(length(idns),1),Metric_base(idns,5),'k');
ids=find(Metric_conflict(:,1)==1);idns=find(Metric_conflict(:,1)==0);
bplot(Metric_conflict(:,5),2,'width',0.5);
scatter(2.*ones(length(ids),1),Metric_conflict(ids,5),'r*');scatter(2.*ones(length(idns),1),Metric_conflict(idns,5),'k');
set(gca,'xtick',1:2,'xticklabel',[{'Baseline'},{'Conflict'}])
ylabel('Neural Features #')


%% RMSE curves
files_base=dir('*Base.mat');

RMSE_base=nan(length(files_base),500);
RMSE_conf=nan(length(files_conf),500);

for ff=1:length(files_base)
    load(files_base(ff).name,'xrmse_curve');
    RMSE_base(ff,1:length(xrmse_curve))=xrmse_curve;
end
ids=find(Metric_base(:,1)==1);idns=find(Metric_base(:,1)==0);
plot(1:500,RMSE_base(idns,:),'k');hold on; 
plot(1:500,nanmean(RMSE_base(idns,:)),'k','linewidth',2);
plot(1:500,RMSE_base(ids,:),'r');
plot(1:500,nanmean(RMSE_base(ids,:)),'r','linewidth',2);
xlim([0 40])
box off

for ff=1:length(files_conf)
    load(files_conf(ff).name,'xrmse_curve');
    RMSE_conf(ff,1:length(xrmse_curve))=xrmse_curve;
end


%% script in order to calculate time frequency analysis of main experiment
% requires Analysis_Corr_SavePreprocessData.m to be run previously



clearvars

%% parameters
F.PathIn                = 'E:\work\data\SSVEP_volmov\EEG\CORR_TFA_RESS_noblinks';
% F.Subjects2Use          = [1:20];
%F.Subjects2Use          = [22 23 25 26 27];
F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
% F.Subjects2Use          = [1 2];


%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    %% preprocessing
    % read in set file
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_corrmat.mat ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    datain = open(sprintf('%s\\VP%02.0f_corrmat.mat',F.PathIn,F.Subjects2Use(i_sub)));
    
    if i_sub == 1
        corr.parameters = datain.corr.parameters;
        corr.parameters_diff = datain.corr.parameters2;
        corr.TFA_baseline = datain.corr.baseline;
        corr.TFA_type = datain.corr.TFAtype;
        corr.outlier_crit = '+-4SD';
        corr.zscdata.xdata=repmat({[]},size(corr.parameters,1),1);
        corr.zscdata.xdata_bc=repmat({[]},size(corr.parameters,1),1);
        corr.zscdata.xdata_diff=repmat({[]},size(corr.parameters_diff,1),1);
        corr.zscdata.ydata=repmat({[]},size(corr.parameters,1),1);
        corr.zscdata.ydata_bc=repmat({[]},size(corr.parameters,1),1);
        corr.zscdata.ydata_diff=repmat({[]},size(corr.parameters_diff,1),1);
        corr.zscdata.subs_xdata=repmat({[]},size(corr.parameters,1),1);
        corr.zscdata.subs_xdata_bc=repmat({[]},size(corr.parameters,1),1);
        corr.zscdata.subs_xdata_diff=repmat({[]},size(corr.parameters_diff,1),1);
        corr.zscdata.subs_ydata=repmat({[]},size(corr.parameters,1),1);
        corr.zscdata.subs_ydata_bc=repmat({[]},size(corr.parameters,1),1);
        corr.zscdata.subs_ydata_diff=repmat({[]},size(corr.parameters_diff,1),1);
        corr.scdata = corr.zscdata;
    end
    
    % extract raw data
    corr.rawdata.xdata(:,i_sub)= datain.corr.xdata;
    corr.rawdata.xdata_bc(:,i_sub)= datain.corr.xdata_bc;
    corr.rawdata.xdata_diff(:,i_sub)= datain.corr.xdata_diff;
    corr.rawdata.ydata(:,i_sub)= datain.corr.ydata;
    corr.rawdata.ydata_bc(:,i_sub)= datain.corr.ydata_bc;
    corr.rawdata.ydata_diff(:,i_sub)= datain.corr.ydata_diff;
    
    % calculate correlations of raw values
    for i_corr = 1:numel(datain.corr.xdata)
        % extract outlier
        t.ind_raw = ...
            datain.corr.xdata{i_corr}>=(mean(datain.corr.xdata{i_corr})+4*std(datain.corr.xdata{i_corr}))|...
            datain.corr.xdata{i_corr}<=(mean(datain.corr.xdata{i_corr})-4*std(datain.corr.xdata{i_corr}))|...
            datain.corr.ydata{i_corr}>=(mean(datain.corr.ydata{i_corr})+4*std(datain.corr.ydata{i_corr}))|...
            datain.corr.ydata{i_corr}<=(mean(datain.corr.ydata{i_corr})-4*std(datain.corr.ydata{i_corr}));
        corr.corr_raw(i_corr).num_outlier_raw(i_sub) = sum(t.ind_raw);
        
        % correlation
        [t.R t.p]=corrcoef(datain.corr.xdata{i_corr}(~t.ind_raw ),datain.corr.ydata{i_corr}(~t.ind_raw ));
        corr.corr_raw(i_corr).R(i_sub) = t.R(2);
        corr.corr_raw(i_corr).p(i_sub) = t.p(2);
        
        % zscore data
        corr.zscdata.xdata{i_corr}=[corr.zscdata.xdata{i_corr}; zscore(datain.corr.xdata{i_corr}(~t.ind_raw ))];
        corr.zscdata.subs_xdata{i_corr}=[corr.zscdata.subs_xdata{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
        corr.zscdata.ydata{i_corr}=[corr.zscdata.ydata{i_corr}; zscore(datain.corr.ydata{i_corr}(~t.ind_raw ))];
        corr.zscdata.subs_ydata{i_corr}=[corr.zscdata.subs_ydata{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
        
        % raw data
        corr.scdata.xdata{i_corr}=[corr.scdata.xdata{i_corr}; datain.corr.xdata{i_corr}(~t.ind_raw )];
        corr.scdata.subs_xdata{i_corr}=[corr.scdata.subs_xdata{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
        corr.scdata.ydata{i_corr}=[corr.scdata.ydata{i_corr}; datain.corr.ydata{i_corr}(~t.ind_raw )];
        corr.scdata.subs_ydata{i_corr}=[corr.scdata.subs_ydata{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
    end
    
    % calculate correlations of bc values
    for i_corr = 1:numel(datain.corr.xdata_bc)
        % extract outlier
        t.ind_raw = ...
            datain.corr.xdata_bc{i_corr}>=(mean(datain.corr.xdata_bc{i_corr})+4*std(datain.corr.xdata_bc{i_corr}))|...
            datain.corr.xdata_bc{i_corr}<=(mean(datain.corr.xdata_bc{i_corr})-4*std(datain.corr.xdata_bc{i_corr}))|...
            datain.corr.ydata_bc{i_corr}>=(mean(datain.corr.ydata_bc{i_corr})+4*std(datain.corr.ydata_bc{i_corr}))|...
            datain.corr.ydata_bc{i_corr}<=(mean(datain.corr.ydata_bc{i_corr})-4*std(datain.corr.ydata_bc{i_corr}));
        corr.corr_bc(i_corr).num_outlier_bc(i_sub) = sum(t.ind_raw);
        
        % correlation
        [t.R t.p]=corrcoef(datain.corr.xdata_bc{i_corr}(~t.ind_raw ),datain.corr.ydata_bc{i_corr}(~t.ind_raw ));
        corr.corr_bc(i_corr).R(i_sub) = t.R(2);
        corr.corr_bc(i_corr).p(i_sub) = t.p(2);
        
        % zscore data
        corr.zscdata.xdata_bc{i_corr}=[corr.zscdata.xdata_bc{i_corr}; zscore(datain.corr.xdata_bc{i_corr}(~t.ind_raw ))];
        corr.zscdata.subs_xdata_bc{i_corr}=[corr.zscdata.subs_xdata_bc{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
        corr.zscdata.ydata_bc{i_corr}=[corr.zscdata.ydata_bc{i_corr}; zscore(datain.corr.ydata_bc{i_corr}(~t.ind_raw ))];
        corr.zscdata.subs_ydata_bc{i_corr}=[corr.zscdata.subs_ydata_bc{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
        
        % raw data
        corr.scdata.xdata_bc{i_corr}=[corr.scdata.xdata_bc{i_corr}; datain.corr.xdata_bc{i_corr}(~t.ind_raw )];
        corr.scdata.subs_xdata_bc{i_corr}=[corr.scdata.subs_xdata_bc{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
        corr.scdata.ydata_bc{i_corr}=[corr.scdata.ydata_bc{i_corr}; datain.corr.ydata_bc{i_corr}(~t.ind_raw )];
        corr.scdata.subs_ydata_bc{i_corr}=[corr.scdata.subs_ydata_bc{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
    end
    
    % calculate correlations of diff
    for i_corr = 1:numel(datain.corr.xdata_diff)
        % extract outlier
        t.ind_raw = ...
            datain.corr.xdata_diff{i_corr}>=(mean(datain.corr.xdata_diff{i_corr})+4*std(datain.corr.xdata_diff{i_corr}))|...
            datain.corr.xdata_diff{i_corr}<=(mean(datain.corr.xdata_diff{i_corr})-4*std(datain.corr.xdata_diff{i_corr}))|...
            datain.corr.ydata_diff{i_corr}>=(mean(datain.corr.ydata_diff{i_corr})+4*std(datain.corr.ydata_diff{i_corr}))|...
            datain.corr.ydata_diff{i_corr}<=(mean(datain.corr.ydata_diff{i_corr})-4*std(datain.corr.ydata_diff{i_corr}));
        corr.corr_diff(i_corr).num_outlier_diff(i_sub) = sum(t.ind_raw);
        
        % correlation
        [t.R t.p]=corrcoef(datain.corr.xdata_diff{i_corr}(~t.ind_raw ),datain.corr.ydata_diff{i_corr}(~t.ind_raw ));
        corr.corr_diff(i_corr).R(i_sub) = t.R(2);
        corr.corr_diff(i_corr).p(i_sub) = t.p(2);
        
        % zscore data
        corr.zscdata.xdata_diff{i_corr}=[corr.zscdata.xdata_diff{i_corr}; zscore(datain.corr.xdata_diff{i_corr}(~t.ind_raw ))];
        corr.zscdata.subs_xdata_diff{i_corr}=[corr.zscdata.subs_xdata_diff{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
        corr.zscdata.ydata_diff{i_corr}=[corr.zscdata.ydata_diff{i_corr}; zscore(datain.corr.ydata_diff{i_corr}(~t.ind_raw ))];
        corr.zscdata.subs_ydata_diff{i_corr}=[corr.zscdata.subs_ydata_diff{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
        
        % raw data
        corr.scdata.xdata_diff{i_corr}=[corr.scdata.xdata_diff{i_corr}; datain.corr.xdata_diff{i_corr}(~t.ind_raw )];
        corr.scdata.subs_xdata_diff{i_corr}=[corr.scdata.subs_xdata_diff{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
        corr.scdata.ydata_diff{i_corr}=[corr.scdata.ydata_diff{i_corr}; datain.corr.ydata_diff{i_corr}(~t.ind_raw )];
        corr.scdata.subs_ydata_diff{i_corr}=[corr.scdata.subs_ydata_diff{i_corr}; repmat(F.Subjects2Use(i_sub),sum(~t.ind_raw),1)];
    end
end

%% plotting
for i_pl = 1:size(corr.parameters,1)
    pl.xdata = corr.zscdata.xdata_bc{i_pl}; pl.labeladd = 'zscore bc';  pl.xlimflag = 1;
    pl.ydata = corr.zscdata.ydata_bc{i_pl};


%     pl.xdata = corr.scdata.xdata_bc{i_pl}; pl.labeladd = 'abs bc';  pl.xlimflag = 1;
%     pl.ydata = corr.scdata.ydata_bc{i_pl};

    
%     pl.xdata = corr.zscdata.xdata{i_pl}; pl.labeladd = 'zscore raw';   pl.xlimflag = 1;
%     pl.ydata = corr.zscdata.ydata{i_pl};

    
%     pl.xdata = corr.scdata.xdata{i_pl}; pl.labeladd = 'abs raw'; pl.xlimflag = 2;
%     pl.ydata = corr.scdata.ydata{i_pl};
    
    figure; scatter(pl.xdata,pl.ydata,'.')
    hold on;
    
    if  pl.xlimflag == 1;
        xlim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
        ylim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
    else
        xlim(1.1*[0 1]*max(abs([pl.xdata;pl.ydata])));
        ylim(1.1*[0 1]*max(abs([pl.xdata;pl.ydata])));
    end
    
    % plot mean
    scatter(mean(pl.xdata),mean(pl.ydata),'co','filled','MarkerEdgeColor','k')
    text(0.9*max(get(gca,'xlim')),mean(pl.ydata),...
            sprintf('[x,y] = [%1.3f,%1.3f]',mean(pl.xdata),mean(pl.ydata)), ...
            'VerticalAlignment','top','HorizontalAlignment','right','FontSize',8,'Color',[1 0 0]);
        
    % correlation
    [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
    if r.P(2)>=.0001
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p = %1.4f\n',r.R(2),r.P(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
    else
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p < .0001\n',r.R(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
    end
    % add fitted line
    r.fit = polyfit(pl.xdata, pl.ydata,1);
    t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
    
    t.pos=get(gca,'Position');
   
    % labels
    xlabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms',pl.labeladd, corr.parameters{i_pl,4},corr.parameters{i_pl,1},corr.parameters{i_pl,2}))
    ylabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms',pl.labeladd, corr.parameters{i_pl,8},corr.parameters{i_pl,5},corr.parameters{i_pl,6}))
     
    grid on
    box on
end

%% density plotting
for i_pl = 1:size(corr.parameters,1)
%     pl.xdata = corr.zscdata.xdata_bc{i_pl}; pl.labeladd = 'zscore bc';  pl.xlimflag = 1;
%     pl.ydata = corr.zscdata.ydata_bc{i_pl};


%     pl.xdata = corr.scdata.xdata_bc{i_pl}; pl.labeladd = 'abs bc';  pl.xlimflag = 1;
%     pl.ydata = corr.scdata.ydata_bc{i_pl};

    
%     pl.xdata = corr.zscdata.xdata{i_pl}; pl.labeladd = 'zscore raw';   pl.xlimflag = 1;
%     pl.ydata = corr.zscdata.ydata{i_pl};

    
    pl.xdata = corr.scdata.xdata{i_pl}; pl.labeladd = 'abs raw'; pl.xlimflag = 2;
    pl.ydata = corr.scdata.ydata{i_pl};
    
    
    figure;
    % define quadratic area
     if  pl.xlimflag == 1;
        t.xlim=[-1 1]*max(abs([pl.xdata;pl.ydata]));
        t.ylim=[-1 1]*max(abs([pl.xdata;pl.ydata]));
    else
        t.xlim=[0 1]*max(abs([pl.xdata;pl.ydata]));
        t.ylim=[0 1]*max(abs([pl.xdata;pl.ydata]));    
    end
    
    
    % run kernel estimation
%     hexscatter(pl.xdata,pl.ydata,'xlim',t.xlim,'ylim',t.ylim,'res',10)
    
    [t.bandwidth,pl.plotmap,t.X,t.Y]=kde2d([pl.xdata,pl.ydata],...
        2^8, [t.xlim(1) t.ylim(1)],[t.xlim(2) t.ylim(2)]);
    imagesc(t.X(1,:),t.Y(:,1),pl.plotmap)
        
%     [pl.plotmap,t.N]=hist3([pl.xdata, pl.ydata],[50 50]);
%     imagesc(t.X(1,:),t.Y(:,1),pl.plotmap)
    
%     dscatter(pl.xdata,pl.ydata,'marker','s','msize',1)
    set(gca,'YDir','normal')
    
    xlim(t.xlim*1/2)
    ylim(t.ylim*1/2)
    
    hold on;
    
    [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
    if r.P(2)>=.0001
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p = %1.4f\n',r.R(2),r.P(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[1 1 1]);
    else
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p < .0001\n',r.R(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[1 1 1]);
    end
    
    % plot mean
    scatter(mean(pl.xdata),mean(pl.ydata),'co','filled','MarkerEdgeColor','k')
    text(0.9*max(max(get(gca,'xlim'))),mean(pl.ydata),...
            sprintf('[x,y] = [%1.3f,%1.3f]',mean(pl.xdata),mean(pl.ydata)), ...
            'VerticalAlignment','top','HorizontalAlignment','right','FontSize',8,'Color',[1 0 0]);
    
    % add fitted line
    r.fit = polyfit(pl.xdata, pl.ydata,1);
    t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
    
    % labels
    xlabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms',pl.labeladd, corr.parameters{i_pl,4},corr.parameters{i_pl,1},corr.parameters{i_pl,2}))
    ylabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms',pl.labeladd, corr.parameters{i_pl,8},corr.parameters{i_pl,5},corr.parameters{i_pl,6}))
      
    colorbar
end

%% plotting diff
for i_pl = 1:size(corr.parameters_diff,1)
    
%     pl.xdata = corr.zscdata.xdata_diff{i_pl}; pl.labeladd = 'zscore';
%     pl.ydata = corr.zscdata.ydata_diff{i_pl};

    
    pl.xdata = corr.scdata.xdata_diff{i_pl}; pl.labeladd = 'abs';
    pl.ydata = corr.scdata.ydata_diff{i_pl};
    
    figure; scatter(pl.xdata,pl.ydata,'.')
    hold on;
    
    xlim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
    ylim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
    
    % plot mean
    scatter(mean(pl.xdata),mean(pl.ydata),'co','filled','MarkerEdgeColor','k')
    text(0.9*max(get(gca,'xlim')),mean(pl.ydata),...
            sprintf('[x,y] = [%1.3f,%1.3f]',mean(pl.xdata),mean(pl.ydata)), ...
            'VerticalAlignment','top','HorizontalAlignment','right','FontSize',8,'Color',[1 0 0]);
        
    % correlation
    [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
    if r.P(2)>=.0001
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p = %1.4f\n',r.R(2),r.P(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
    else
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p < .0001\n',r.R(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
    end
    % add fitted line
    r.fit = polyfit(pl.xdata, pl.ydata,1);
    t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
    
    t.pos=get(gca,'Position');
        
    xlabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms - [%1.0f %1.0f]ms',...
        pl.labeladd, corr.parameters_diff{i_pl,5},corr.parameters_diff{i_pl,1},corr.parameters_diff{i_pl,2},corr.parameters_diff{i_pl,3}))
    ylabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms - [%1.0f %1.0f]ms',...
        pl.labeladd, corr.parameters_diff{i_pl,10},corr.parameters_diff{i_pl,6},corr.parameters_diff{i_pl,7}, corr.parameters_diff{i_pl,8}))
    
    grid on
    box on
end

%% density plotting diff
for i_pl = 1:size(corr.parameters_diff,1)
%     pl.xdata = corr.zscdata.xdata_diff{i_pl}; pl.labeladd = 'zscore';
%     pl.ydata = corr.zscdata.ydata_diff{i_pl};

    
    pl.xdata = corr.scdata.xdata_diff{i_pl}; pl.labeladd = 'abs';
    pl.ydata = corr.scdata.ydata_diff{i_pl};
    
    
    figure;
    % define quadratic area
    t.xlim=[-1 1]*max(abs([pl.xdata;pl.ydata]));
    t.ylim=[-1 1]*max(abs([pl.xdata;pl.ydata]));
    
    
    % run kernel estimation
%     hexscatter(pl.xdata,pl.ydata,'xlim',t.xlim,'ylim',t.ylim,'res',10)
    
    [t.bandwidth,pl.plotmap,t.X,t.Y]=kde2d([pl.xdata,pl.ydata],...
        2^8, [t.xlim(1) t.ylim(1)],[t.xlim(2) t.ylim(2)]);
    imagesc(t.X(1,:),t.Y(:,1),pl.plotmap)
        
%     [pl.plotmap,t.N]=hist3([pl.xdata, pl.ydata],[50 50]);
%     imagesc(t.X(1,:),t.Y(:,1),pl.plotmap)
    
%     dscatter(pl.xdata,pl.ydata,'marker','s','msize',1)
    set(gca,'YDir','normal')
    
    xlim(t.xlim*1/2)
    ylim(t.ylim*1/2)
    
    hold on;
    
    [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
    if r.P(2)>=.0001
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p = %1.4f\n',r.R(2),r.P(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[1 1 1]);
    else
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p < .0001\n',r.R(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[1 1 1]);
    end
    
    % plot mean
    scatter(mean(pl.xdata),mean(pl.ydata),'co','filled','MarkerEdgeColor','k')
    text(0.9*max(max(get(gca,'xlim'))),mean(pl.ydata),...
            sprintf('[x,y] = [%1.3f,%1.3f]',mean(pl.xdata),mean(pl.ydata)), ...
            'VerticalAlignment','top','HorizontalAlignment','right','FontSize',8,'Color',[1 0 0]);
    
    % add fitted line
    r.fit = polyfit(pl.xdata, pl.ydata,1);
    t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
    
    % labels
   xlabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms - [%1.0f %1.0f]ms',...
        pl.labeladd, corr.parameters_diff{i_pl,5},corr.parameters_diff{i_pl,1},corr.parameters_diff{i_pl,2},corr.parameters_diff{i_pl,3}))
    ylabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms - [%1.0f %1.0f]ms',...
        pl.labeladd, corr.parameters_diff{i_pl,10},corr.parameters_diff{i_pl,6},corr.parameters_diff{i_pl,7}, corr.parameters_diff{i_pl,8}))  
    colorbar
end


%% plotting (experimental: regress linear relationship)
for i_pl = 1:size(corr.parameters,1)
%     pl.xdata = corr.zscdata.xdata_bc{i_pl}; pl.labeladd = 'zscore bc';  pl.xlimflag = 1;
%     r.fit = polyfit(corr.scdata.xdata{i_pl}, corr.scdata.ydata{i_pl},1);
%     pl.ydata = corr.zscdata.ydata_bc{i_pl}-(r.fit(1).*corr.zscdata.xdata_bc{i_pl});


%     pl.xdata = corr.scdata.xdata_bc{i_pl}; pl.labeladd = 'abs bc';  pl.xlimflag = 1;
%     r.fit = polyfit(corr.scdata.xdata{i_pl}, corr.scdata.ydata{i_pl},1);
%     pl.ydata = corr.scdata.ydata_bc{i_pl}-(r.fit(1).*corr.scdata.xdata_bc{i_pl});

    
%     pl.xdata = corr.zscdata.xdata{i_pl}; pl.labeladd = 'zscore raw';   pl.xlimflag = 1;
%     r.fit = polyfit(corr.scdata.xdata{i_pl}, corr.scdata.ydata{i_pl},1);
%     pl.ydata = corr.zscdata.ydata{i_pl}-(r.fit(1).*corr.zscdata.xdata{i_pl});

    
    pl.xdata = corr.scdata.xdata{i_pl}; pl.labeladd = 'abs raw'; pl.xlimflag = 2;
    r.fit = polyfit(corr.scdata.xdata{i_pl}, corr.scdata.ydata{i_pl},1);
    pl.ydata = corr.scdata.ydata{i_pl}-(r.fit(1).*corr.scdata.xdata{i_pl});
    
    figure; scatter(pl.xdata,pl.ydata,'.')
    hold on;
    
    if  pl.xlimflag == 1;
        xlim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
        ylim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
    else
        xlim(1.1*[0 1]*max(abs([pl.xdata;pl.ydata])));
        ylim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
    end
    
    % plot mean
    scatter(mean(pl.xdata),mean(pl.ydata),'co','filled','MarkerEdgeColor','k')
    text(0.9*max(get(gca,'xlim')),mean(pl.ydata),...
            sprintf('[x,y] = [%1.3f,%1.3f]',mean(pl.xdata),mean(pl.ydata)), ...
            'VerticalAlignment','top','HorizontalAlignment','right','FontSize',8,'Color',[1 0 0]);
        
    % correlation
    [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
    if r.P(2)>=.0001
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p = %1.4f\n',r.R(2),r.P(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
    else
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p < .0001\n',r.R(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
    end
    % add fitted line
    r.fit = polyfit(pl.xdata, pl.ydata,1);
    t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
    
    t.pos=get(gca,'Position');
   
    % labels
    xlabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms',pl.labeladd, corr.parameters{i_pl,4},corr.parameters{i_pl,1},corr.parameters{i_pl,2}))
    ylabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms',pl.labeladd, corr.parameters{i_pl,8},corr.parameters{i_pl,5},corr.parameters{i_pl,6}))
     
    grid on
    box on
    
    % add counts
    text(0.8*min(get(gca,'xlim')),0.9*max(get(gca,'ylim')),...
        sprintf('N = %1.0f of %1.0f\np = %1.4f',...
        sum((pl.xdata<0)&(pl.ydata>0)), numel(pl.xdata), sum(pl.xdata<0&pl.ydata>0)/numel(pl.xdata)), ...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[0 0 0]);
    text(0.8*max(get(gca,'xlim')),0.9*max(get(gca,'ylim')),...
        sprintf('N = %1.0f of %1.0f\np = %1.4f',...
        sum((pl.xdata>0)&(pl.ydata>0)), numel(pl.xdata), sum(pl.xdata>0&pl.ydata>0)/numel(pl.xdata)), ...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[0 0 0]);
    text(0.8*min(get(gca,'xlim')),0.8*min(get(gca,'ylim')),...
        sprintf('N = %1.0f of %1.0f\np = %1.4f',...
        sum((pl.xdata<0)&(pl.ydata<0)), numel(pl.xdata), sum(pl.xdata<0&pl.ydata<0)/numel(pl.xdata)), ...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[0 0 0]);
    text(0.8*max(get(gca,'xlim')),0.8*min(get(gca,'ylim')),...
        sprintf('N = %1.0f of %1.0f\np = %1.4f',...
        sum((pl.xdata>0)&(pl.ydata<0)), numel(pl.xdata), sum(pl.xdata>0&pl.ydata<0)/numel(pl.xdata)), ...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[0 0 0]);
end

%% experimental plotting (experimental: regress linear relationship) | density with histograms
pl.vernum = 5; pl.hnum = 4;
pl.sbpln_hist1 = [1 2 3];
pl.sbpln_hist2 = [8 12 16];
pl.sbpln_main = [5 6 7 9 10 11 13 14 15];
pl.sbpln_cbar = [17 18 19];
for i_pl = 1:size(corr.parameters,1)
    pl.xdata = corr.zscdata.xdata_bc{i_pl}; pl.labeladd = 'zscore bc';  pl.xlimflag = 1;
    r.fit = polyfit(corr.scdata.xdata{i_pl}, corr.scdata.ydata{i_pl},1);
    pl.ydata = corr.zscdata.ydata_bc{i_pl}-(r.fit(1).*corr.zscdata.xdata_bc{i_pl});


%     pl.xdata = corr.scdata.xdata_bc{i_pl}; pl.labeladd = 'abs bc';  pl.xlimflag = 1;
%     r.fit = polyfit(corr.scdata.xdata{i_pl}, corr.scdata.ydata{i_pl},1);
%     pl.ydata = corr.scdata.ydata_bc{i_pl}-(r.fit(1).*corr.scdata.xdata_bc{i_pl});

    
%     pl.xdata = corr.zscdata.xdata{i_pl}; pl.labeladd = 'zscore raw';   pl.xlimflag = 1;
%     r.fit = polyfit(corr.scdata.xdata{i_pl}, corr.scdata.ydata{i_pl},1);
%     pl.ydata = corr.zscdata.ydata{i_pl}-(r.fit(1).*corr.zscdata.xdata{i_pl});

    
%     pl.xdata = corr.scdata.xdata{i_pl}; pl.labeladd = 'abs raw'; pl.xlimflag = 2;
%     r.fit = polyfit(corr.scdata.xdata{i_pl}, corr.scdata.ydata{i_pl},1);
%     pl.ydata = corr.scdata.ydata{i_pl}-(r.fit(1).*corr.scdata.xdata{i_pl});
    
    figure;
    h.sp(i_pl,1)=subplot(pl.vernum,pl.hnum,pl.sbpln_main);
    % define quadratic area
    t.xlim=[-1 1]*max(abs([pl.xdata;pl.ydata]));
    t.ylim=[-1 1]*max(abs([pl.xdata;pl.ydata]));
    
    
    % run kernel estimation
%     hexscatter(pl.xdata,pl.ydata,'xlim',t.xlim,'ylim',t.ylim,'res',10)
    
    [t.bandwidth,pl.plotmap,t.X,t.Y]=kde2d([pl.xdata,pl.ydata],...
        2^8, [t.xlim(1) t.ylim(1)],[t.xlim(2) t.ylim(2)]);
    imagesc(t.X(1,:),t.Y(:,1),pl.plotmap)
        
%     [pl.plotmap,t.N]=hist3([pl.xdata, pl.ydata],[50 50]);
%     imagesc(t.X(1,:),t.Y(:,1),pl.plotmap)
    
%     dscatter(pl.xdata,pl.ydata,'marker','s','msize',1)
    set(gca,'YDir','normal')
    
    xlim(t.xlim*1/2)
    ylim(t.ylim*1/2)
    
    hold on;
    
    [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
    if r.P(2)>=.0001
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p = %1.4f\n',r.R(2),r.P(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[1 1 1]);
    else
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p < .0001\n',r.R(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[1 1 1]);
    end
    
    % plot mean
    scatter(mean(pl.xdata),mean(pl.ydata),'co','filled','MarkerEdgeColor','k')
    text(0.9*max(max(get(gca,'xlim'))),mean(pl.ydata),...
            sprintf('[x,y] = [%1.3f,%1.3f]',mean(pl.xdata),mean(pl.ydata)), ...
            'VerticalAlignment','top','HorizontalAlignment','right','FontSize',8,'Color',[1 0 0]);
    
    % add fitted line
    r.fit = polyfit(pl.xdata, pl.ydata,1);
    t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
    
    % labels
    xlabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms',pl.labeladd, corr.parameters{i_pl,4},corr.parameters{i_pl,1},corr.parameters{i_pl,2}))
    ylabel(sprintf('%s %s | [%1.0f %1.0f]Hz | [%1.0f %1.0f]ms',pl.labeladd, corr.parameters{i_pl,8},corr.parameters{i_pl,5},corr.parameters{i_pl,6}))
    
    grid on
    box on
     
    % histograms
    h.sp(i_pl,2)=subplot(pl.vernum,pl.hnum,pl.sbpln_hist1)
    [pl.counts,pl.bins] = hist(pl.xdata,linspace(t.xlim(1)*1/2,t.xlim(2)*1/2,60));
    bar(pl.bins,pl.counts)
    xlim(t.xlim*1/2)
    vline(mean(pl.xdata),'c')
    
    h.sp(i_pl,3)=subplot(pl.vernum,pl.hnum,pl.sbpln_hist2)
    [pl.counts,pl.bins] = hist(pl.ydata,linspace(t.xlim(1)*1/2,t.xlim(2)*1/2,60));
    barh(pl.bins,pl.counts)
    ylim(t.ylim*1/2)
    hline(mean(pl.ydata),'c')
    
    %
    h.sp(i_pl,4)=subplot(pl.vernum,pl.hnum,pl.sbpln_cbar)
    caxis([0 max(max(pl.plotmap))])
    colorbar('Location','southoutside')
    set(h.sp(i_pl,4),'Visible','off');
end

%% experimental plotting diff
pl.vernum = 5; pl.hnum = 4;
pl.sbpln_hist1 = [1 2 3];
pl.sbpln_hist2 = [8 12 16];
pl.sbpln_main = [5 6 7 9 10 11 13 14 15];
pl.sbpln_cbar = [17 18 19];
for i_pl = 1:size(corr.parameters_diff,1)
%     pl.xdata = corr.zscdata.xdata_diff{i_pl}; pl.labeladd = 'zscore';
%     pl.ydata = corr.zscdata.ydata_diff{i_pl};
    
    pl.xdata = corr.scdata.xdata_diff{i_pl}; pl.labeladd = 'abs';
    pl.ydata = corr.scdata.ydata_diff{i_pl};
    
    figure;
    h.sp(i_pl,1)=subplot(pl.vernum,pl.hnum,pl.sbpln_main);
    % define quadratic area
    t.xlim=[-1 1]*max(abs([pl.xdata;pl.ydata]));
    t.ylim=[-1 1]*max(abs([pl.xdata;pl.ydata]));
    
    
    % run kernel estimation
%     hexscatter(pl.xdata,pl.ydata,'xlim',t.xlim,'ylim',t.ylim,'res',10)
    
    [t.bandwidth,pl.plotmap,t.X,t.Y]=kde2d([pl.xdata,pl.ydata],...
        2^8, [t.xlim(1) t.ylim(1)],[t.xlim(2) t.ylim(2)]);
    imagesc(t.X(1,:),t.Y(:,1),pl.plotmap)
        
%     [pl.plotmap,t.N]=hist3([pl.xdata, pl.ydata],[50 50]);
%     imagesc(t.X(1,:),t.Y(:,1),pl.plotmap)
    
%     dscatter(pl.xdata,pl.ydata,'marker','s','msize',1)
    set(gca,'YDir','normal')
    
    xlim(t.xlim*1/2)
    ylim(t.ylim*1/2)
    
    hold on;
    
    [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
    if r.P(2)>=.0001
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p = %1.4f\n',r.R(2),r.P(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[1 1 1]);
    else
        text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
            sprintf('R = %1.3f p < .0001\n',r.R(2)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8,'Color',[1 1 1]);
    end
    
    % plot mean
    scatter(mean(pl.xdata),mean(pl.ydata),'co','filled','MarkerEdgeColor','k')
    text(0.9*max(max(get(gca,'xlim'))),mean(pl.ydata),...
            sprintf('[x,y] = [%1.3f,%1.3f]',mean(pl.xdata),mean(pl.ydata)), ...
            'VerticalAlignment','top','HorizontalAlignment','right','FontSize',8,'Color',[1 0 0]);
    
    % add fitted line
    r.fit = polyfit(pl.xdata, pl.ydata,1);
    t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
    
    % labels
    xlabel(sprintf('%s %s | [%1.0f %1.0f]Hz\n[%1.0f %1.0f]ms - [%1.0f %1.0f]ms',...
        pl.labeladd, corr.parameters_diff{i_pl,5},corr.parameters_diff{i_pl,1},corr.parameters_diff{i_pl,2},corr.parameters_diff{i_pl,3}))
    ylabel(sprintf('%s %s | [%1.0f %1.0f]Hz\n[%1.0f %1.0f]ms - [%1.0f %1.0f]ms',...
        pl.labeladd, corr.parameters_diff{i_pl,10},corr.parameters_diff{i_pl,6},corr.parameters_diff{i_pl,7}, corr.parameters_diff{i_pl,8}))
    
    % histograms
    h.sp(i_pl,2)=subplot(pl.vernum,pl.hnum,pl.sbpln_hist1)
    [pl.counts,pl.bins] = hist(pl.xdata,linspace(t.xlim(1)*1/2,t.xlim(2)*1/2,60));
    bar(pl.bins,pl.counts)
    xlim(t.xlim*1/2)
    vline(mean(pl.xdata),'c')
    
    h.sp(i_pl,3)=subplot(pl.vernum,pl.hnum,pl.sbpln_hist2)
    [pl.counts,pl.bins] = hist(pl.ydata,linspace(t.xlim(1)*1/2,t.xlim(2)*1/2,60));
    barh(pl.bins,pl.counts)
    ylim(t.ylim*1/2)
    hline(mean(pl.ydata),'c')
    
    %
    h.sp(i_pl,4)=subplot(pl.vernum,pl.hnum,pl.sbpln_cbar)
    caxis([0 max(max(pl.plotmap))])
    colorbar('Location','southoutside')
    set(h.sp(i_pl,4),'Visible','off');
end



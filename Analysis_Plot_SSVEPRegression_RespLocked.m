%% script to analyze
% based on RESS components



clearvars

%% parameters
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_lagged_regression_resplocked';
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_lagged_regression_resplocked_alpha';

F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 25 26 27];
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 17 18 19 20 25 26 27];


pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% load data

% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_lagreg_withSSVEP_noblinks ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    DATAIn=load(sprintf('%s\\VP%02.0f_exp_lagreg_withSSVEP_noblinks.mat',F.PathIn,F.Subjects2Use(i_sub)));
    
    % preallocate memory
    if i_sub == 1
        REG.data_b = nan([size(DATAIn.REG.data_b) numel(numel(F.Subjects2Use))]);
        REG.data_p = nan([size(DATAIn.REG.data_p) numel(numel(F.Subjects2Use))]);
        REG.data_tstat = nan([size(DATAIn.REG.data_tstat) numel(numel(F.Subjects2Use))]);
        REG.data_Rsqr_adj = nan([size(DATAIn.REG.data_Rsqr_adj) numel(numel(F.Subjects2Use))]);
        REG.lagtime = DATAIn.REG.time;
        REG.resptime = DATAIn.REG.resptime;
        REG.freq = DATAIn.REG.freq;
        REG.chanlocs = DATAIn.REG.chanlocs;
        REG.params = DATAIn.REG.params;
    end
    
    % extract data
    REG.params(i_sub)=DATAIn.REG.params;
    REG.data_b(:,:,:,:,i_sub)=DATAIn.REG.data_b;
    REG.data_p(:,:,:,:,i_sub)=DATAIn.REG.data_p;
    REG.data_tstat(:,:,:,:,i_sub)=DATAIn.REG.data_tstat;
    REG.data_Rsqr_adj(:,:,:,:,i_sub)=DATAIn.REG.data_Rsqr_adj;
    
%     % extract number of trials post-hoc
%     t.tail2P = nan(1,100000);
%     for i_n=50000:200000
%         tail2P(i_n) = 2*tcdf(-abs(DATAIn.REG.data_tstat(100)),(i_n));
%     end
%     [find(tail2P>DATAIn.REG.data_p(100),1,'first') find(tail2P<DATAIn.REG.data_p(100),1,'last')-1]
%     x = -5:0.01:5;
%     y1 = tpdf(x,5);
%     y2 = tpdf(x,15);
%     y3 = tpdf(x,5000);
%     z = normpdf(x,0,1);
%     plot(x,y1,'-.',x,y2,'--',x,y3,'.',x,z,'-')
%     legend('Student''s t Distribution with \nu=5', ...
%         'Student''s t Distribution with \nu=25', ...
%         'Student''s t Distribution with \nu=5000', ...
%         'Standard Normal Distribution','Location','best')
%     title('Student''s t and Standard Normal pdfs')
    
    
end

%% plot data

%% plot data for separate signals | different estimation of significance!

% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [8 12]; pl.ROI = 'visual alpha';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'visual SSVEP';
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [10 14]; pl.ROI = 'motor alpha';
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [18 30]; pl.ROI = 'motor beta';
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [16 23]; pl.ROI = 'parietal beta';
pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [8 12]; pl.ROI = 'parietal alpha';
% pl.elec2plot = {'FP1';'FP2'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'frontal SSVEP';

pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({REG.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freq2plot_i = dsearchn(REG.freq', pl.freq2plot');

pl.plotdata = []; pl.clims = []; pl.titadd = {};
pl.plotdata(:,:,1) = squeeze(mean(mean(mean(REG.data_b(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2),5)); 
pl.plotdata(:,:,2) = squeeze(mean(mean(mean(REG.data_tstat(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2),5));
pl.plotdata(:,:,3) = squeeze(mean(mean(mean(REG.data_Rsqr_adj(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2),5));

% average of p-values
pl.plotdata(:,:,4) = squeeze(abs(log10(mean(mean(mean(REG.data_p(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2),5))));

% directly test t-values to be smaller than icdf('Normal',0.05,0,1) and larger than  icdf('Normal',0.95,0,1)
[tt.h_l,tt.p_l,tt.ci_l,tt.stats_l] = ttest(squeeze(mean(mean(REG.data_tstat(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2)), ...
    icdf('Normal',0.025,0,1),'Dim',3,'Alpha', 0.05/2,'Tail','left');
[tt.h_r,tt.p_r,tt.ci_r,tt.stats_r] = ttest(squeeze(mean(mean(REG.data_tstat(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2)), ...
    icdf('Normal',0.975,0,1),'Dim',3,'Alpha', 0.05/2,'Tail','right');
tt.idx = tt.p_l<tt.p_r;
tt.p_lr = abs(log10(tt.p_l)); tt.p_lr(~tt.idx)=abs(log10(tt.p_r(~tt.idx)));
pl.plotdata(:,:,5) = tt.p_lr;
pl.titadd{5} ='log p-val | ttest on tvals (|tval| >= 1.96)';

% directly test beta values
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(mean(mean(REG.data_b(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2)),0,'Dim',3);
pl.plotdata(:,:,6) = abs(log10(tt.p));
pl.titadd{6} ='log p-val | ttest on betas (|beta| >= 0)';

% plot only certain time to responses or time lags
pl.lag2plot = REG.lagtime([1 end]);
% pl.lag2plot = REG.lagtime([4 end-4]);
pl.lag2plot_i = dsearchn(REG.lagtime',pl.lag2plot');

% pl.respt2plot = REG.resptime([1 end]);
pl.respt2plot = REG.resptime([3 end-3]);
pl.respt2plot_i = dsearchn(REG.resptime',pl.respt2plot');

pl.plotdata = pl.plotdata(pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:);


% color limits
pl.clims(1,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,1)))); pl.titadd{1} ='betas';
pl.clims(2,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,2)))); pl.titadd{2} ='t-stat';
pl.clims(3,:) = [0 1]*max(max(abs(pl.plotdata(:,:,3)))); pl.titadd{3} ='r^2 adj';
pl.tdata = pl.plotdata(:,:,4); pl.clims(4,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata))))); pl.titadd{4} ='log p-val | mean(p)';
pl.tdata = pl.plotdata(:,:,5); pl.clims(5,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata)))));
pl.tdata = pl.plotdata(:,:,6); pl.clims(6,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata)))));


figure;
for i_d = 1:6
    subplot(2,3,i_d)
    imagesc(REG.resptime(pl.respt2plot_i(1):pl.respt2plot_i(2)),REG.lagtime(pl.lag2plot_i(1):pl.lag2plot_i(2)),...
        pl.plotdata(:,:,i_d),pl.clims(i_d,:))
    set(gca,'YDir','normal')
    title(sprintf('%s\ntime lagged regression on %s signal', pl.titadd{i_d}, REG.params(1).TFA_regressor_signal), 'FontSize',8,'Interpreter','none')
    h.cb2=colorbar;
    xlabel('time to response')
    ylabel('lag in ms')
    set(gca,'FontSize',8)
    if i_d == 1 | i_d == 2
        colormap(gca,flipud(cbrewer2('RdBu')))
    elseif i_d >= 4 % for pvalues
        t.pcriterion = abs(log10(0.05));
        if max(max(pl.plotdata(:,:,i_d)))<t.pcriterion
            % temp.colormap = repmat([0.5 0.5 0.5],100,1);
            t.colormap = repmat(linspace(1,0.3,1000)',1,3);
        else
            t.border = ceil((t.pcriterion / max(max(pl.plotdata(:,:,i_d))))*1000);
            % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
            t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
        end
        colormap(gca,t.colormap) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
        t.yticks = get(h.cb2,'YTick');
        set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<pl.clims(i_d,end),1,'last')), ...
            'YTickLabel',pl.plegend(1:find(pl.pcorrect<pl.clims(i_d,end),1,'last')),'FontSize',8)
        
    elseif i_d == 3
        colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    end
end
h.ax2 = axes('Position', [-0.02, .44, .15, .15],'Visible','off');
% topoplot(find(pl.elec2plot_i),REG.chanlocs(1:64),'style','blank','electrodes', 'on','whitebk','on');
topoplot(find(pl.elec2plot_i),REG.chanlocs(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});
title(pl.ROI)

%% plot data for separate signals | different estimation of significance! | flipped

% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [8 12]; pl.ROI = 'visual alpha';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'visual SSVEP';
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [10 14]; pl.ROI = 'motor alpha';
pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [18 30]; pl.ROI = 'motor beta';
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [16 23]; pl.ROI = 'parietal beta';
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [8 12]; pl.ROI = 'parietal alpha';
% pl.elec2plot = {'FP1';'FP2'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'frontal SSVEP';

pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({REG.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freq2plot_i = dsearchn(REG.freq', pl.freq2plot');

pl.plotdata = []; pl.clims = []; pl.titadd = {};
pl.plotdata(:,:,1) = squeeze(mean(mean(mean(REG.data_b(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2),5))'; 
pl.plotdata(:,:,2) = squeeze(mean(mean(mean(REG.data_tstat(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2),5))';
pl.plotdata(:,:,3) = squeeze(mean(mean(mean(REG.data_Rsqr_adj(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2),5))';

% average of p-values
pl.plotdata(:,:,4) = squeeze(abs(log10(mean(mean(mean(REG.data_p(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2),5))))';

% directly test t-values to be smaller than icdf('Normal',0.05,0,1) and larger than  icdf('Normal',0.95,0,1)
[tt.h_l,tt.p_l,tt.ci_l,tt.stats_l] = ttest(squeeze(mean(mean(REG.data_tstat(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2)), ...
    icdf('Normal',0.025,0,1),'Dim',3,'Alpha', 0.05/2,'Tail','left');
[tt.h_r,tt.p_r,tt.ci_r,tt.stats_r] = ttest(squeeze(mean(mean(REG.data_tstat(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2)), ...
    icdf('Normal',0.975,0,1),'Dim',3,'Alpha', 0.05/2,'Tail','right');
tt.idx = tt.p_l<tt.p_r;
tt.p_lr = abs(log10(tt.p_l)); tt.p_lr(~tt.idx)=abs(log10(tt.p_r(~tt.idx)));
pl.plotdata(:,:,5) = tt.p_lr';
pl.titadd{5} ='log p-val | ttest on tvals (|tval| >= 1.96)';

% directly test beta values
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(mean(mean(REG.data_b(pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),:,:,:),1),2)),0,'Dim',3);
pl.plotdata(:,:,6) = abs(log10(tt.p))';
pl.titadd{6} ='log p-val | ttest on betas (|beta| >= 0)';

% plot only certain time to responses or time lags
pl.lag2plot = REG.lagtime([1 end]);
% pl.lag2plot = REG.lagtime([4 end-4]);
pl.lag2plot_i = dsearchn(REG.lagtime',pl.lag2plot');

% pl.respt2plot = REG.resptime([1 end]);
pl.respt2plot = REG.resptime([3 end-3]);
pl.respt2plot_i = dsearchn(REG.resptime',pl.respt2plot');

pl.plotdata = pl.plotdata(pl.respt2plot_i(1):pl.respt2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),:);


% color limits
pl.clims(1,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,1)))); pl.titadd{1} ='betas';
pl.clims(2,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,2)))); pl.titadd{2} ='t-stat';
pl.clims(3,:) = [0 1]*max(max(abs(pl.plotdata(:,:,3)))); pl.titadd{3} ='r^2 adj';
pl.tdata = pl.plotdata(:,:,4); pl.clims(4,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata))))); pl.titadd{4} ='log p-val | mean(p)';
pl.tdata = pl.plotdata(:,:,5); pl.clims(5,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata)))));
pl.tdata = pl.plotdata(:,:,6); pl.clims(6,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata)))));


figure;
for i_d = 1:6
    subplot(2,3,i_d)
    imagesc(REG.lagtime(pl.lag2plot_i(1):pl.lag2plot_i(2)),REG.resptime(pl.respt2plot_i(1):pl.respt2plot_i(2)),...
        pl.plotdata(:,:,i_d),pl.clims(i_d,:))
    set(gca,'YDir','normal')
    title(sprintf('%s\ntime lagged regression on %s signal', pl.titadd{i_d}, REG.params(1).TFA_regressor_signal), 'FontSize',8,'Interpreter','none')
    h.cb2=colorbar;
    ylabel('time to response')
    xlabel('lag in ms')
    set(gca,'FontSize',8)
    if i_d == 1 | i_d == 2
        colormap(gca,flipud(cbrewer2('RdBu')))
    elseif i_d >= 4 % for pvalues
        t.pcriterion = abs(log10(0.05));
        if max(max(pl.plotdata(:,:,i_d)))<t.pcriterion
            % temp.colormap = repmat([0.5 0.5 0.5],100,1);
            t.colormap = repmat(linspace(1,0.3,1000)',1,3);
        else
            t.border = ceil((t.pcriterion / max(max(pl.plotdata(:,:,i_d))))*1000);
            % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
            t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
        end
        colormap(gca,t.colormap) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
        t.yticks = get(h.cb2,'YTick');
        set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<pl.clims(i_d,end),1,'last')), ...
            'YTickLabel',pl.plegend(1:find(pl.pcorrect<pl.clims(i_d,end),1,'last')),'FontSize',8)
        
    elseif i_d == 3
        colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    end
end
h.ax2 = axes('Position', [-0.02, .44, .15, .15],'Visible','off');
% topoplot(find(pl.elec2plot_i),REG.chanlocs(1:64),'style','blank','electrodes', 'on','whitebk','on');
topoplot(find(pl.elec2plot_i),REG.chanlocs(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});
title(pl.ROI)




%% plot averages across all lags and across all response times
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [8 12]; pl.ROI = 'visual alpha';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'visual SSVEP';
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [10 14]; pl.ROI = 'motor alpha';
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [18 30]; pl.ROI = 'motor beta';
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [16 23]; pl.ROI = 'parietal beta';
pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [8 12]; pl.ROI = 'parietal alpha';
% pl.elec2plot = {'FP1';'FP2'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'frontal SSVEP';

pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({REG.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% frequency to plot
pl.freq2plot_i = dsearchn(REG.freq', pl.freq2plot');

% plot only certain time to responses or time lags
pl.lag2plot = REG.lagtime([1 end]);
% pl.lag2plot = REG.lagtime([4 end-4]);
pl.lag2plot_i = dsearchn(REG.lagtime',pl.lag2plot');

% pl.respt2plot = REG.resptime([1 end]);
pl.respt2plot = REG.resptime([3 end-3]);
pl.respt2plot_i = dsearchn(REG.resptime',pl.respt2plot');

% across lags
pl.plotdata_lag = []; pl.titadd = {};
pl.plotdata_lag(:,1) = squeeze(mean(mean(mean(mean(REG.data_b(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),3),5));
pl.plotdata_lag(:,2) = squeeze(mean(mean(mean(mean(REG.data_tstat(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),3),5));
pl.plotdata_lag(:,3) = squeeze(mean(mean(mean(mean(REG.data_Rsqr_adj(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),3),5));

% directly test beta values
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(mean(mean(mean(REG.data_b(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),3)),0,'Dim',2);
% pl.plotdata_lag(:,4) = abs(log10(tt.p));
pl.plotdata_lag(:,4) = tt.p;


% across response times
pl.plotdata_resp = [];
pl.plotdata_resp(:,1) = squeeze(mean(mean(mean(mean(REG.data_b(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),4),5)); 
pl.plotdata_resp(:,2) = squeeze(mean(mean(mean(mean(REG.data_tstat(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),4),5));
pl.plotdata_resp(:,3) = squeeze(mean(mean(mean(mean(REG.data_Rsqr_adj(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),4),5));

% directly test beta values
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(mean(mean(mean(REG.data_b(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),4)),0,'Dim',2);
% pl.plotdata_resp(:,4) = abs(log10(tt.p));
pl.plotdata_resp(:,4) = tt.p;

% title add
% color limits
pl.titadd = {'betas','t-stat','r^2 adj','log p-val | ttest on betas (|beta| >= 0)'};
pl.titadd = {'betas','t-stat','r^2 adj','p-val | ttest on betas (|beta| >= 0)'};

% plot data
fig = figure;
set(gcf,'Position',[100 100 900 900],'PaperPositionMode','auto')

for i_dat = 1:4
    subplot(4,2,1+2*(i_dat-1))
    
    % plot data averaged across lags
    plot(REG.resptime(pl.respt2plot_i(1):pl.respt2plot_i(2)),pl.plotdata_lag(:,i_dat))
    xlim(REG.resptime([pl.respt2plot_i(1) pl.respt2plot_i(2)]))
    % title
    title(sprintf('M_l_a_g: %s',pl.titadd{i_dat}))
    grid on
%     if i_dat == 4
%         set(gca,'yscale','log','ydir','reverse')
%         set(gca,'YTickLabel',flipud(get(gca,'YTick')))
% %         set(gca,'yscale','log','ydir','reverse','YTickLabel',get(gca,'YTick'))
%         xlabel('time to response in ms')
%     end
    if i_dat == 4
        set(gca,'yscale','log','ydir','reverse')
        ylim([min(get(gca,'YLim')) 1])
        set(gca,'YTickLabel',flipud(get(gca,'YTick')))
%         set(gca,'yscale','log','ydir','reverse','YTickLabel',get(gca,'YTick'))
        xlabel('time to response in ms')
        
        if min(get(gca,'YLim'))
            hold on
            plot(get(gca,'XLim'),[.05 .05],'b')
        end
    end
    set(gca, 'Position',get(gca, 'Position')+[0.05 0 0 0])
    
    
    subplot(4,2,2+2*(i_dat-1))
    
    % plot data averaged across lags
    plot(REG.lagtime(pl.lag2plot_i(1):pl.lag2plot_i(2)),pl.plotdata_resp(:,i_dat))
    xlim(REG.lagtime([pl.lag2plot_i(1) pl.lag2plot_i(2)]))
    % title
    title(sprintf('M_r_e_s_p: %s',pl.titadd{i_dat}))
    grid on
    if i_dat == 4
        set(gca,'yscale','log','ydir','reverse')
        ylim([min(get(gca,'YLim')) 1])
        set(gca,'YTickLabel',flipud(get(gca,'YTick')))
%         set(gca,'yscale','log','ydir','reverse','YTickLabel',get(gca,'YTick'))
        xlabel('lag time in ms')
        
        if min(get(gca,'YLim'))
            hold on
            plot(get(gca,'XLim'),[.05 .05],'b')
        end
    end
    set(gca, 'Position',get(gca, 'Position')+[0.05 0 0 0])

end

h.a1 = axes('position',[0.025 0.45 0.12 0.12],'Visible','off');
topoplot(find(pl.elec2plot_i),REG.chanlocs(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});
title(sprintf('%s\n[%1.1f %1.1f] Hz',pl.ROI,pl.freq2plot),'FontSize',8)


%% plot averages across defined lags
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [8 12]; pl.ROI = 'visual alpha';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'visual SSVEP';
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [10 14]; pl.ROI = 'motor alpha';
pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [18 30]; pl.ROI = 'motor beta';
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [16 23]; pl.ROI = 'parietal beta';
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [8 12]; pl.ROI = 'parietal alpha';
% pl.elec2plot = {'FP1';'FP2'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'frontal SSVEP';

pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({REG.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% frequency to plot
pl.freq2plot_i = dsearchn(REG.freq', pl.freq2plot');

% plot only certain time to responses or time lags
pl.lag2plot = [0 500]
pl.lag2plot = [-200 0];;
% pl.lag2plot = REG.lagtime([4 end-4]);
pl.lag2plot_i = dsearchn(REG.lagtime',pl.lag2plot');

% pl.respt2plot = REG.resptime([1 end]);
pl.respt2plot = REG.resptime([3 end-3]);
pl.respt2plot_i = dsearchn(REG.resptime',pl.respt2plot');

% across lags
pl.plotdata_lag = []; pl.titadd = {};
pl.plotdata_lag(:,1) = squeeze(mean(mean(mean(mean(REG.data_b(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),3),5));
pl.plotdata_lag(:,2) = squeeze(mean(mean(mean(mean(REG.data_tstat(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),3),5));
pl.plotdata_lag(:,3) = squeeze(mean(mean(mean(mean(REG.data_Rsqr_adj(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),3),5));

% directly test beta values
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(mean(mean(mean(REG.data_b(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),3)),0,'Dim',2);
% pl.plotdata_lag(:,4) = abs(log10(tt.p));
pl.plotdata_lag(:,4) = tt.p;


% title add
pl.titadd = {'betas','t-stat','r^2 adj','p-val | ttest on betas (|beta| >= 0)'};

% plot data
fig = figure;
set(gcf,'Position',[100 100 700 900],'PaperPositionMode','auto')

for i_dat = 1:4
    subplot(4,1,i_dat)
    
    % plot data averaged across lags
    plot(REG.resptime(pl.respt2plot_i(1):pl.respt2plot_i(2)),pl.plotdata_lag(:,i_dat))
    xlim(REG.resptime([pl.respt2plot_i(1) pl.respt2plot_i(2)]))
    % title
    title(sprintf('%s',pl.titadd{i_dat}))
    grid on
%     if i_dat == 4
%         set(gca,'yscale','log','ydir','reverse')
%         set(gca,'YTickLabel',flipud(get(gca,'YTick')))
% %         set(gca,'yscale','log','ydir','reverse','YTickLabel',get(gca,'YTick'))
%         xlabel('time to response in ms')
%     end
    if i_dat == 4
        set(gca,'yscale','log','ydir','reverse')
        ylim([min(get(gca,'YLim')) 1])
        set(gca,'YTickLabel',flipud(get(gca,'YTick')))
%         set(gca,'yscale','log','ydir','reverse','YTickLabel',get(gca,'YTick'))
        xlabel('time to response in ms')
        
        if min(get(gca,'YLim'))
            hold on
            plot(get(gca,'XLim'),[.05 .05],'b')
        end
    end
    set(gca, 'Position',get(gca, 'Position')+[0.07 0 -0.02 0])

end

h.a1 = axes('position',[0.025 0.45 0.12 0.17],'Visible','off');
topoplot(find(pl.elec2plot_i),REG.chanlocs(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});
title(sprintf('%s\n[%1.1f %1.1f] Hz\nfor lag [%1.0f %1.0f] ms',pl.ROI,pl.freq2plot,pl.lag2plot),'FontSize',8)


%% plot averages across defined times to response
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [8 12]; pl.ROI = 'visual alpha';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'visual SSVEP';
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [10 14]; pl.ROI = 'motor alpha';
pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.freq2plot = [18 30]; pl.ROI = 'motor beta';
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [16 23]; pl.ROI = 'parietal beta';
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.freq2plot = [8 12]; pl.ROI = 'parietal alpha';
% pl.elec2plot = {'FP1';'FP2'}; pl.freq2plot = [14.167 14.167]; pl.ROI = 'frontal SSVEP';

pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({REG.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% frequency to plot
pl.freq2plot_i = dsearchn(REG.freq', pl.freq2plot');

% plot only certain time to responses or time lags
pl.lag2plot = REG.lagtime([1 end]);
% pl.lag2plot = REG.lagtime([4 end-4]);
pl.lag2plot_i = dsearchn(REG.lagtime',pl.lag2plot');

pl.respt2plot = [0 1000];
% pl.respt2plot = [-1000 0];
% pl.respt2plot = [-1800 1800];
pl.respt2plot_i = dsearchn(REG.resptime',pl.respt2plot');

% across response times
pl.plotdata_resp = [];
pl.plotdata_resp(:,1) = squeeze(mean(mean(mean(mean(REG.data_b(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),4),5)); 
pl.plotdata_resp(:,2) = squeeze(mean(mean(mean(mean(REG.data_tstat(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),4),5));
pl.plotdata_resp(:,3) = squeeze(mean(mean(mean(mean(REG.data_Rsqr_adj(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),4),5));

% directly test beta values
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(mean(mean(mean(REG.data_b(...
    pl.elec2plot_i,pl.freq2plot_i(1):pl.freq2plot_i(2),pl.lag2plot_i(1):pl.lag2plot_i(2),pl.respt2plot_i(1):pl.respt2plot_i(2),:),1),2),4)),0,'Dim',2);
% pl.plotdata_resp(:,4) = abs(log10(tt.p));
pl.plotdata_resp(:,4) = tt.p;

% title add
pl.titadd = {'betas','t-stat','r^2 adj','p-val | ttest on betas (|beta| >= 0)'};

% plot data
fig = figure;
set(gcf,'Position',[100 100 600 900],'PaperPositionMode','auto')

for i_dat = 1:4
    subplot(4,1,i_dat)
    
    % plot data averaged across lags
    plot(REG.lagtime(pl.lag2plot_i(1):pl.lag2plot_i(2)),pl.plotdata_resp(:,i_dat))
    xlim(REG.lagtime([pl.lag2plot_i(1) pl.lag2plot_i(2)]))
    % title
    title(sprintf('M_r_e_s_p: %s',pl.titadd{i_dat}))
    grid on
    if i_dat == 4
        set(gca,'yscale','log','ydir','reverse')
        ylim([min(get(gca,'YLim')) 1])
        set(gca,'YTickLabel',flipud(get(gca,'YTick')))
%         set(gca,'yscale','log','ydir','reverse','YTickLabel',get(gca,'YTick'))
        xlabel('lag time in ms')
        
        if min(get(gca,'YLim'))
            hold on
            plot(get(gca,'XLim'),[.05 .05],'b')
        end
    end
    set(gca, 'Position',get(gca, 'Position')+[0.07 0 -0.02 0])

end

h.a1 = axes('position',[0.025 0.45 0.12 0.17],'Visible','off');
topoplot(find(pl.elec2plot_i),REG.chanlocs(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});
title(sprintf('%s\n[%1.1f %1.1f] Hz',pl.ROI,pl.freq2plot),'FontSize',8)
title(sprintf('%s\n[%1.1f %1.1f] Hz\nfor resp.time [%1.0f %1.0f] ms',pl.ROI,pl.freq2plot,pl.respt2plot),'FontSize',8)
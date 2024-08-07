%% script in order to calculate regression based on SSVEPs
% based on RESS components



clearvars

%% parameters
F.PathIn                = 'E:\work\data\SSVEP_volmov\EEG\TFA_lagged_regression';

F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 25 26 27];
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 17 18 19 20 25 26 27];


pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% load data

% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
%     fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_lagreg_withSSVEP ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
%     DATAIn=load(sprintf('%s\\VP%02.0f_exp_lagreg_withSSVEP',F.PathIn,F.Subjects2Use(i_sub)));
%     fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_lagreg_withSSVEP_nonlinks ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
%     DATAIn=load(sprintf('%s\\VP%02.0f_exp_lagreg_withSSVEP_noblinks.mat',F.PathIn,F.Subjects2Use(i_sub)));
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_lagreg_withSSVEP_nonlinks_v2 ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    DATAIn=load(sprintf('%s\\VP%02.0f_exp_lagreg_withSSVEP_noblinks_v2.mat',F.PathIn,F.Subjects2Use(i_sub)));
%     fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_lagreg_withSSVEP_nonlinks_v2_detrended ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
%     DATAIn=load(sprintf('%s\\VP%02.0f_exp_lagreg_withSSVEP_noblinks_v2_detrended.mat',F.PathIn,F.Subjects2Use(i_sub)));
    
    % preallocate memory
    if i_sub == 1
        REG.data_b = nan([size(DATAIn.REG.data_b) numel(numel(F.Subjects2Use))]);
        REG.data_p = nan([size(DATAIn.REG.data_p) numel(numel(F.Subjects2Use))]);
        REG.data_tstat = nan([size(DATAIn.REG.data_tstat) numel(numel(F.Subjects2Use))]);
        REG.data_Rsqr_adj = nan([size(DATAIn.REG.data_Rsqr_adj) numel(numel(F.Subjects2Use))]);
        REG.time = DATAIn.REG.time;
        REG.freq = DATAIn.REG.freq;
        REG.chanlocs = DATAIn.REG.chanlocs;
        REG.params = DATAIn.REG.params;
    end
    
    % extract data
    REG.params(i_sub)=DATAIn.REG.params;
    REG.data_b(:,:,:,i_sub)=DATAIn.REG.data_b;
    REG.data_p(:,:,:,i_sub)=DATAIn.REG.data_p;
    REG.data_tstat(:,:,:,i_sub)=DATAIn.REG.data_tstat;
    REG.data_Rsqr_adj(:,:,:,i_sub)=DATAIn.REG.data_Rsqr_adj;
    
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
%% single plots
% version 1
% pl.elec2plot = {'Oz'};
pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; pl.ROI = 'steady state'; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; pl.ROI = 'steady state'; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; pl.ROI = 'vis alpha';% vis alpha
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.ROI = 'vis alpha'; % vis alpha II
% pl.elec2plot = {'POz'}; pl.ROI = 'vis alpha'; % vis alpha II
% pl.elec2plot = {'POz';'Oz'}; pl.ROI = 'vis alpha'; % vis alpha III
% pl.elec2plot = {'C3';'CP3'}; pl.ROI = 'motor'; % motor alpha/beta
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.ROI = 'motor'; % motor alpha/beta II
% pl.elec2plot = {'F1';'F2';'Fz'}; pl.ROI = 'frontal'; % frontal
% pl.elec2plot = {'FP1';'FP2'}; pl.ROI = 'frontal'; % frontal
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.ROI = 'parietal'; % parietal
% pl.elec2plot = {'F5';'F3';'F1';'FC3'}; pl.ROI = 'pre-motor'; % pre-motor
% pl.elec2plot = {'FC3';'FC3';'FC1';'FCz';'Fz';'F5';'F3';'F1';'AF7';'AF3';'AFz';'FP1';'FPz'}; pl.ROI = 'left-centro-frontal'; % pre-motor
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({REG.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freq2plot = REG.freq([1 end]);
% pl.freq2plot = [REG.freq(1) 12];
pl.f_idx = dsearchn(REG.freq', pl.freq2plot');

pl.plotdata = []; pl.clims = []; pl.titadd = {};
pl.plotdata(:,:,1) = squeeze(mean(mean(REG.data_b(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1),4)); 
pl.clims(1,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,1)))); pl.titadd{1} ='betas';
pl.plotdata(:,:,2) = squeeze(mean(mean(REG.data_tstat(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1),4));
pl.clims(2,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,2)))); pl.titadd{2} ='t-stat';

% average of p-values
% pl.plotdata(:,:,3) = squeeze(abs(log10(mean(mean(REG.data_p(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1),4)))); pl.titadd{3} ='log p-val';

% directly test t-values to be smaller than icdf('Normal',0.05,0,1) and larger than  icdf('Normal',0.95,0,1)
[tt.h_l,tt.p_l,tt.ci_l,tt.stats_l] = ttest(squeeze(mean(REG.data_tstat(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1)), icdf('Normal',0.025,0,1),...
    'Dim',3,'Alpha', 0.05/2,'Tail','left');
[tt.h_r,tt.p_r,tt.ci_r,tt.stats_r] = ttest(squeeze(mean(REG.data_tstat(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1)), icdf('Normal',0.975,0,1),...
    'Dim',3,'Alpha', 0.05/2,'Tail','right');
tt.idx = tt.p_l<tt.p_r;
tt.p_lr = abs(log10(tt.p_l)); tt.p_lr(~tt.idx)=abs(log10(tt.p_r(~tt.idx)));
pl.plotdata(:,:,3) = tt.p_lr;
pl.titadd{3} ='log p-val';
pl.tdata = pl.plotdata(:,:,3); pl.clims(3,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata)))));

pl.plotdata(:,:,4) = squeeze(mean(mean(REG.data_Rsqr_adj(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1),4)); 
pl.clims(4,:) = [0 1]*max(max(abs(pl.plotdata(:,:,4)))); pl.titadd{4} ='r^2 adj';

figure;
for i_d = 1:4
    subplot(2,2,i_d)
    imagesc(REG.time,REG.freq(pl.f_idx(1):pl.f_idx(2)),pl.plotdata(:,:,i_d),pl.clims(i_d,:))
    set(gca,'YDir','normal')
    title(sprintf('%s\ntime lagged regression on SSVEP signal | %s', pl.titadd{i_d}, pl.ROI), 'FontSize',8)
    h.cb2=colorbar;
    xlabel('lag in ms')
    ylabel('frequency in Hz')
    hline(14.16667,'m')
    set(gca,'FontSize',8)
    if i_d == 1 | i_d == 2
        colormap(gca,flipud(cbrewer2('RdBu')))
    elseif i_d == 3 % for pvalues
        t.pcriterion = abs(log10(0.05));
        if max(pl.tdata(:))<t.pcriterion
            % temp.colormap = repmat([0.5 0.5 0.5],100,1);
            t.colormap = repmat(linspace(1,0.3,1000)',1,3);
        else
            t.border = ceil((t.pcriterion / max(pl.tdata(:)))*1000);
            % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
            t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
        end
        colormap(gca,t.colormap) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
        t.yticks = get(h.cb2,'YTick');
        set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<pl.clims(i_d,end),1,'last')), ...
            'YTickLabel',pl.plegend(1:find(pl.pcorrect<pl.clims(i_d,end),1,'last')),'FontSize',8)
        
    elseif i_d == 4
        colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    end
end
h.ax2 = axes('Position', [0.43, .44, .15, .15],'Visible','off');
topoplot(find(pl.elec2plot_i),REG.chanlocs(1:64),'style','blank','electrodes', 'on','whitebk','on');


%% single plots | different estimation of significance!
% version 1
% pl.elec2plot = {'Oz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; pl.ROI = 'steady state'; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; pl.ROI = 'steady state'; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; pl.ROI = 'vis alpha';% vis alpha
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.ROI = 'vis alpha'; % vis alpha II
% pl.elec2plot = {'POz'}; pl.ROI = 'vis alpha'; % vis alpha II
% pl.elec2plot = {'POz';'Oz'}; pl.ROI = 'vis alpha'; % vis alpha III
% pl.elec2plot = {'C3';'CP3'}; pl.ROI = 'motor'; % motor alpha/beta
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; pl.ROI = 'motor'; % motor alpha/beta II
% pl.elec2plot = {'F1';'F2';'Fz'}; pl.ROI = 'frontal'; % frontal
% pl.elec2plot = {'FP1';'FP2'}; pl.ROI = 'frontal'; % frontal
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; pl.ROI = 'parietal'; % parietal
% pl.elec2plot = {'F5';'F3';'F1';'FC3'}; pl.ROI = 'pre-motor'; % pre-motor
% pl.elec2plot = {'FC3';'FC3';'FC1';'FCz';'Fz';'F5';'F3';'F1';'AF7';'AF3';'AFz';'FP1';'FPz'}; pl.ROI = 'left-centro-frontal'; % pre-motor
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({REG.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% pl.freq2plot = REG.freq([1 end]);
pl.freq2plot = [REG.freq(1) 31];
pl.freq2plot = [3.7 REG.freq(end)];
% pl.freq2plot = [3.7 31];
pl.f_idx = dsearchn(REG.freq', pl.freq2plot');

pl.plotdata = []; pl.clims = []; pl.titadd = {};
pl.plotdata(:,:,1) = squeeze(mean(mean(REG.data_b(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1),4)); 
pl.clims(1,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,1)))); pl.titadd{1} ='betas';
pl.plotdata(:,:,2) = squeeze(mean(mean(REG.data_tstat(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1),4));
pl.clims(2,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,2)))); pl.titadd{2} ='t-stat';
pl.plotdata(:,:,3) = squeeze(mean(mean(REG.data_Rsqr_adj(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1),4)); 
pl.clims(3,:) = [0 1]*max(max(abs(pl.plotdata(:,:,3)))); pl.titadd{3} ='r^2 adj';

% average of p-values
pl.plotdata(:,:,4) = squeeze(abs(log10(mean(mean(REG.data_p(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1),4)))); 
pl.tdata = pl.plotdata(:,:,4); pl.clims(4,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata))))); pl.titadd{4} ='log p-val | mean(p)';

% directly test t-values to be smaller than icdf('Normal',0.05,0,1) and larger than  icdf('Normal',0.95,0,1)
[tt.h_l,tt.p_l,tt.ci_l,tt.stats_l] = ttest(squeeze(mean(REG.data_tstat(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1)), icdf('Normal',0.025,0,1),...
    'Dim',3,'Alpha', 0.05/2,'Tail','left');
[tt.h_r,tt.p_r,tt.ci_r,tt.stats_r] = ttest(squeeze(mean(REG.data_tstat(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1)), icdf('Normal',0.975,0,1),...
    'Dim',3,'Alpha', 0.05/2,'Tail','right');
tt.idx = tt.p_l<tt.p_r;
tt.p_lr = abs(log10(tt.p_l)); tt.p_lr(~tt.idx)=abs(log10(tt.p_r(~tt.idx)));
pl.plotdata(:,:,5) = tt.p_lr;
pl.titadd{5} ='log p-val | ttest on tvals (|tval| >= 1.96)';
pl.tdata = pl.plotdata(:,:,5); pl.clims(5,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata)))));

% directly test beta values
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(mean(REG.data_b(pl.elec2plot_i,pl.f_idx(1):pl.f_idx(2),:,:),1)),0,'Dim',3);
pl.plotdata(:,:,6) = abs(log10(tt.p));
pl.titadd{6} ='log p-val | ttest on betas (|beta| >= 0)';
pl.tdata = pl.plotdata(:,:,6); pl.clims(6,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata)))));



figure;
for i_d = 1:6
    subplot(2,3,i_d)
    imagesc(REG.time,REG.freq(pl.f_idx(1):pl.f_idx(2)),pl.plotdata(:,:,i_d),pl.clims(i_d,:))
    set(gca,'YDir','normal')
    title(sprintf('%s\ntime lagged regression on SSVEP signal | %s', pl.titadd{i_d}, pl.ROI), 'FontSize',8)
    h.cb2=colorbar;
    xlabel('lag in ms')
    ylabel('frequency in Hz')
    hline(14.16667,'m')
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
topoplot(find(pl.elec2plot_i),REG.chanlocs(1:64),'style','blank','electrodes', 'on','whitebk','on');


%% topographical plot
pl.plotdata = []; pl.clims = []; pl.titadd = {};
pl.plotdata(:,:,:,1) = mean(REG.data_b,4); pl.clims(1,:) = [-1 1]*max(max(max(abs(pl.plotdata(:,:,:,1))))); pl.titadd{1} ='betas';
pl.plotdata(:,:,:,2) = mean(REG.data_tstat,4); pl.clims(2,:) = [-1 1]*max(max(max(abs(pl.plotdata(:,:,:,2))))); pl.titadd{2} ='t-stat';

% averaging p-values?
% pl.plotdata(:,:,:,3) = abs(log10(mean(REG.data_p,4))); pl.titadd{3} ='log p-val';
% pl.tdata = pl.plotdata(:,:,:,3); pl.clims(3,:) = [0 1]*max(max(max(abs(pl.tdata(~isinf(pl.tdata))))));

% testing t-values against 1.96
% [tt.h_l,tt.p_l,tt.ci_l,tt.stats_l] = ttest(squeeze(REG.data_tstat), icdf('Normal',0.025,0,1),...
%     'Dim',4,'Alpha', 0.05/2,'Tail','left');
% [tt.h_r,tt.p_r,tt.ci_r,tt.stats_r] = ttest(squeeze(REG.data_tstat), icdf('Normal',0.975,0,1),...
%     'Dim',4,'Alpha', 0.05/2,'Tail','right');
% tt.idx = tt.p_l<tt.p_r;
% tt.p_lr = abs(log10(tt.p_l)); tt.p_lr(~tt.idx)=abs(log10(tt.p_r(~tt.idx)));
% pl.plotdata(:,:,:,3) = tt.p_lr;
% pl.titadd{3} ='log p-val';
% pl.tdata = pl.plotdata(:,:,:,3); pl.clims(3,:) = [0 1]*max(max(max(abs(pl.tdata(~isinf(pl.tdata))))));

% testing betas against 0

% directly test beta values
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(REG.data_b),0,'Dim',4);
pl.plotdata(:,:,:,3) = abs(log10(tt.p));
pl.titadd{6} ='log p-val | (|beta| >= 0)';
l.tdata = pl.plotdata(:,:,:,3); pl.clims(3,:) = [0 1]*max(max(max(abs(pl.tdata(~isinf(pl.tdata))))));



pl.plotdata(:,:,:,4) = mean(REG.data_Rsqr_adj,4); pl.titadd{4} ='r^2 adj';
pl.clims(4,:) = [0 0.005]; %pl.clims(4,:) = [0 1]*max(max(max(abs(pl.plotdata(:,:,:,4)))));
pl.flims = [REG.freq(1) 31];
% pl.flims = [REG.freq(1) REG.freq(end)];


for i_d = 1:4
    if i_d == 3
        plot_data_topoarray(REG.chanlocs, pl.plotdata(:,:,:,i_d),'TFA','times',REG.time,...
            'freqs',REG.freq,'title', pl.titadd{i_d},'clim', pl.clims(i_d,:),'ylim_f',pl.flims,'p_vals',true)
    else
        plot_data_topoarray(REG.chanlocs, pl.plotdata(:,:,:,i_d),'TFA','times',REG.time,...
            'freqs',REG.freq,'title', pl.titadd{i_d},'clim', pl.clims(i_d,:),'ylim_f',pl.flims)
    end
end


%% plot data course
%% single plots
% version 1
pl.elec2plot = {'O1';'Oz';'O2';'POz'}; pl.ROI = 'visual central'; pl.freq2plot = [10 14];
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.ROI = 'vis alpha'; pl.freq2plot = [8 12];
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.ROI = 'vis Theta'; pl.freq2plot = [4 7];
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.ROI = 'vis SSVEP'; pl.freq2plot = [14.167 14.167];
% pl.elec2plot = {'C3';'CP3'}; pl.freq2plot = [10 14]; pl.ROI = 'mot alpha';
% pl.elec2plot = {'FP1';'FPz';'FP2'}; pl.ROI = 'frontal eyes';
% pl.elec2plot = {'F1';'Fz';'F2';'FC1';'FCz';'FC2'}; pl.freq2plot = [8 12]; pl.ROI = 'frontal alpha';
% pl.elec2plot = {'F1';'Fz';'F2';'FC1';'FCz';'FC2'}; pl.freq2plot = [3 7]; pl.ROI = 'frontal theta';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({REG.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqidx = dsearchn(REG.freq',pl.freq2plot');
pl.plotdata = []; pl.titadd = {};
pl.plotdata(:,:,1) = squeeze(mean(mean(REG.data_b(pl.elec2plot_i,pl.freqidx(1):pl.freqidx(2),:,:),1),2)); 
pl.titadd{1} ='betas';
pl.plotdata(:,:,2) = squeeze(mean(mean(REG.data_tstat(pl.elec2plot_i,pl.freqidx(1):pl.freqidx(2),:,:),1),2));
pl.titadd{2} ='t-stat';

% summing p-values
% pl.plotdata(:,:,3) = squeeze(abs(log10(mean(mean((REG.data_p(pl.elec2plot_i,pl.freqidx(1):pl.freqidx(2),:,:)),1),2))));
% pl.titadd{3} ='log p-val';

% running t-tests on t-stat values
[tt.h_l,tt.p_l,tt.ci_l,tt.stats_l] = ttest(squeeze(REG.data_tstat(:,pl.f_idx(1):pl.f_idx(2),:,:)), icdf('Normal',0.025,0,1),...
    'Dim',4,'Alpha', 0.05/2,'Tail','left');
[tt.h_r,tt.p_r,tt.ci_r,tt.stats_r] = ttest(squeeze(REG.data_tstat(:,pl.f_idx(1):pl.f_idx(2),:,:)), icdf('Normal',0.975,0,1),...
    'Dim',4,'Alpha', 0.05/2,'Tail','right');
tt.idx = tt.p_l<tt.p_r;
tt.p_lr = abs(log10(tt.p_l)); tt.p_lr(~tt.idx)=abs(log10(tt.p_r(~tt.idx)));
pl.plotdata(:,:,:,3) = tt.p_lr;
pl.titadd{3} ='log p-val';

pl.plotdata(:,:,4) = squeeze(mean(mean(REG.data_Rsqr_adj(pl.elec2plot_i,pl.freqidx(1):pl.freqidx(2),:,:),1),2)); 
pl.titadd{4} ='r^2 adj';

figure;
subplot(4,1,1)
plot(REG.time,mean(pl.plotdata(:,:,1),2));
ylabel(pl.titadd{1});
set(gca,'XTickLabel',[])
title(sprintf('%s | %s | [%1.1f %1.1f]Hz',pl.ROI,vararg2str(pl.elec2plot),pl.freq2plot))
hline(0,'k')
subplot(4,1,2)
plot(REG.time,mean(pl.plotdata(:,:,2),2));
ylabel(pl.titadd{2});
set(gca,'XTickLabel',[])
hline(0,'k')
subplot(4,1,3)
plot(REG.time,mean(pl.plotdata(:,:,3),2));
ylabel(pl.titadd{3});
set(gca,'XTickLabel',[])
hline(abs(log10(.05)),'r'); hline(abs(log10(.01)),'b'); hline(abs(log10(.001)),'g')
subplot(4,1,4)
plot(REG.time,mean(pl.plotdata(:,:,4),2));
ylabel(pl.titadd{4});

%% plot previously calculated TFA data
% different RESS settings for VP 17 and VP 19

clearvars
%% parameters
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\EOG_Data\';
% F.Subjects2Use          = [1:4];
F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27];

% TFA.baseline            = [-2750 -2250];
% TFA.baseline            = [-2250 -2000];
EOG.baseline            = [-3500 -3000];
EOG.baseline            = [-4000 -3500];
EOG.baseline            = [-4500 -4000];



pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% loop for subjects
for i_sub=1:numel(F.Subjects2Use)
    %%
    % read in TFA data
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_BLNKRESP.mat ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    
    temp.eog = open(sprintf('%sVP%02.0f_exp_BLNKRESP.mat',F.PathIn,F.Subjects2Use(i_sub)));
    
    % a little troubleshooting
    temp.eog.BLNK.df2 = cell2struct([...
        num2cell(temp.eog.BLNK.df.onset_samples); num2cell(temp.eog.BLNK.df.onset_times);...
        num2cell(temp.eog.BLNK.df.onset_samples_intrial); num2cell(temp.eog.BLNK.df.onset_times_intrial);...
        num2cell(temp.eog.BLNK.df.onset_trialindex); num2cell(temp.eog.BLNK.df.onset_type);...
        num2cell(repmat(F.Subjects2Use(i_sub),1,numel(temp.eog.BLNK.df.onset_samples)))],...
        [fieldnames(temp.eog.BLNK.df);{'sub'}]);
        
    % preallocate memory
    if i_sub == 1
        DATA_df =  temp.eog.BLNK.df2
    else
        DATA_df = [DATA_df;  temp.eog.BLNK.df2];
    end
    
    % extract data
    temp.datfields = fieldnames(temp.eog.BLNK);
    temp.datfield_idx = [1:18 20];
    for i_dat = 1:numel(temp.datfield_idx)
        DATA(i_sub).(temp.datfields{temp.datfield_idx(i_dat)})=...
            temp.eog.BLNK.(temp.datfields{temp.datfield_idx(i_dat)});
    end
    
    clear temp
    
end

%% do permutation test of uniform distribution of blinks and responses
% and visualize data
perm.timewin = [0,min([DATA.totaltime])/1000]; % time window for which blinks are anaysed
perm.bins = 60;
perm.n_reps = 10000; % number of repetitions
clear t
t.blink_idx = strcmp([DATA_df.onset_type],'blink');
t.resp_idx = strcmp([DATA_df.onset_type],'resp');

% first extract data across subjects
for i_sub = 1:numel(F.Subjects2Use)
    % only blinks in range are relevant
    t.blnk_data = ...
        [DATA_df(([DATA_df.onset_times]/1000)<perm.timewin(2) & [DATA_df.sub]==F.Subjects2Use(i_sub) & t.blink_idx).onset_times]/1000;
    t.resp_data = ...
        [DATA_df(([DATA_df.onset_times]/1000)<perm.timewin(2) & [DATA_df.sub]==F.Subjects2Use(i_sub) & t.resp_idx).onset_times]/1000;
    
    [t.blnk_hist.N,t.blnk_hist.edges,t.blnk_hist.bin] = histcounts(t.blnk_data,perm.bins,'BinLimits',perm.timewin);
    [t.resp_hist.N,t.resp_hist.edges,t.resp_hist.bin] = histcounts(t.resp_data,perm.bins,'BinLimits',perm.timewin);
    
    perm.bincenter = t.resp_hist.edges(1:end-1)+(diff(t.resp_hist.edges)./2);
    perm.o_data_blnk(i_sub,:)=t.blnk_hist.N;
    perm.o_data_resp(i_sub,:)=t.resp_hist.N;
end
% check plotting
% figure; bar(perm.bincenter,sum(perm.o_data_blnk))
% figure; bar(perm.bincenter,sum(perm.o_data_resp))

% do permutation
perm.sh_data_blnk = nan(perm.n_reps,perm.bins);
perm.sh_data_resp = nan(perm.n_reps,perm.bins);
for i_rep = 1:perm.n_reps
    % shuffle index for all subjects
    t.idx = cell2mat(cellfun(@(x) randperm(x), repmat(num2cell(perm.bins),1,numel(F.Subjects2Use)),'UniformOutput',false)');
    % create histogram according to this shuffling
    perm.sh_data_blnk(i_rep,:) = sum(perm.o_data_blnk(t.idx));
    perm.sh_data_resp(i_rep,:) = sum(perm.o_data_resp(t.idx));
end

% find median and percentiles for each bin of empirical distribution
perm.sh_percentiles_blnk = nan(3,perm.bins);
perm.sh_percentiles_resp = nan(3,perm.bins);
for i_bin = 1:perm.bins
%     perm.sh_percentiles_blnk(:,i_bin) = prctile(perm.sh_data_blnk(:,i_bin), [2.5 50 97.5]);
%     perm.sh_percentiles_resp(:,i_bin) = prctile(perm.sh_data_resp(:,i_bin), [2.5 50 97.5]);
    perm.sh_percentiles_blnk(:,i_bin) = prctile(perm.sh_data_blnk(:,i_bin), [2.5/perm.bins 50 100-2.5/perm.bins]);
    perm.sh_percentiles_resp(:,i_bin) = prctile(perm.sh_data_resp(:,i_bin), [2.5/perm.bins 50 100-2.5/perm.bins]);
end

% visualize data

%%%%%%%%%%%%
% manually [blinks and responses together]
figs{1} = figure;
set(gcf,'Position',[100 100 600 200],'PaperPositionMode','auto')
subplot(4,1,1)
t.respidx = strcmp([DATA_df.onset_type],'resp');
h.l1 = line([[DATA_df(t.respidx).onset_times]./1000; [DATA_df(t.respidx).onset_times]./1000],...
    [ones(1,sum(t.respidx))+0.025; ones(1,sum(t.respidx))+0.5], ...
    'Color', [0.8 0.1 0.1 0.05],'LineWidth',1.5);
hold on
t.blnkidx = strcmp([DATA_df.onset_type],'blink');
h.l2 = line([[DATA_df(t.blnkidx).onset_times]./1000; [DATA_df(t.blnkidx).onset_times]./1000],...
    [ones(1,sum(t.blnkidx))-0.5; ones(1,sum(t.blnkidx))-0.025], ...
    'Color', [0.2 0.4 1 0.05],'LineWidth',1.5);
set(gca,'YTickLabels',[],'YTick',[],'XTickLabels',[],'XTick',[],'Box','on')
xlim([0,min([DATA.totaltime])/1000])
ylim([0.45,1.55])
% plot vertical lines
h.l3 = line(repmat(linspace(1/6,5/6,5).*min([DATA.totaltime])/1000,2,1),repmat(get(gca,'YLim'),5,1)',...
    'Color', [0.2 0.6 0.2 1],'LineWidth',1.5);
%     xlabel('time in s')
%     legend([h.l1(1) h.l2(1)], {'finger taps';'blinks'},'Location','southoutside','Orientation','horizontal')
title('distribution of events across experimental run')

subplot(4,1,[2 3 4])
h.hist1 = histogram([DATA_df(t.respidx).onset_times]./1000,60,'BinLimits',[0,min([DATA.totaltime])/1000],...
    'FaceColor', [0.8 0.1 0.1]);
hold on
h.hist2 = histogram([DATA_df(t.blnkidx).onset_times]./1000,60,'BinLimits',[0,min([DATA.totaltime])/1000],...
    'FaceColor', [0.2 0.4 1]);
h.li = plot(perm.bincenter,[perm.sh_percentiles_blnk; perm.sh_percentiles_resp]);
set(h.li(1),'Color', [0.2 0.4 1], 'LineWidth',1.5, 'LineStyle',':')
set(h.li(2),'Color', [0.2 0.4 1], 'LineWidth',1.5, 'LineStyle','-')
set(h.li(3),'Color', [0.2 0.4 1], 'LineWidth',1.5, 'LineStyle',':')
set(h.li(4),'Color', [0.8 0.1 0.1], 'LineWidth',1.5, 'LineStyle',':')
set(h.li(5),'Color', [0.8 0.1 0.1], 'LineWidth',1.5, 'LineStyle','-')
set(h.li(6),'Color', [0.8 0.1 0.1], 'LineWidth',1.5, 'LineStyle',':')
%    title('distribution of events across experimental run')
xlim([0,min([DATA.totaltime])/1000])
ylim(get(gca,'YLim'))
h.l4 = line(repmat(linspace(1/6,5/6,5).*min([DATA.totaltime])/1000,2,1),repmat(get(gca,'YLim'),5,1)',...
    'Color', [0.2 0.6 0.2 1],'LineWidth',1.5);
legend([h.hist1 h.hist2],{'finger taps';'blinks'},'Location','southoutside','Orientation','horizontal')
xlabel('time in s')

%%%%%%%%%%%%%%%
% [blinks and responses separately]
figs{2} = figure;
set(gcf,'Position',[100 100 600 200],'PaperPositionMode','auto')
subplot(4,1,1)
t.respidx = strcmp([DATA_df.onset_type],'resp');
h.l1 = line([[DATA_df(t.respidx).onset_times]./1000; [DATA_df(t.respidx).onset_times]./1000],...
    [ones(1,sum(t.respidx)); ones(1,sum(t.respidx))+0.5], ...
    'Color', [0.8 0.1 0.1 0.05],'LineWidth',1.5);
hold on
set(gca,'YTickLabels',[],'YTick',[],'XTickLabels',[],'XTick',[],'Box','on')
xlim([0,min([DATA.totaltime])/1000])
ylim([0.95,1.55])
% plot vertical lines
h.l3 = line(repmat(linspace(1/6,5/6,5).*min([DATA.totaltime])/1000,2,1),repmat(get(gca,'YLim'),5,1)',...
    'Color', [0.2 0.6 0.2 1],'LineWidth',1.5);
%     xlabel('time in s')
%     legend([h.l1(1) h.l2(1)], {'finger taps';'blinks'},'Location','southoutside','Orientation','horizontal')
title('distribution of responses across experimental run')

subplot(4,1,[2 3 4])
h.hist1 = histogram([DATA_df(t.respidx).onset_times]./1000,60,'BinLimits',[0,min([DATA.totaltime])/1000],...
    'FaceColor', [0.8 0.1 0.1]);
hold on
h.li = plot(perm.bincenter,[perm.sh_percentiles_resp]);
set(h.li(1),'Color', [0.8 0.1 0.1], 'LineWidth',1.5, 'LineStyle',':')
set(h.li(2),'Color', [0.8 0.1 0.1], 'LineWidth',1.5, 'LineStyle','-')
set(h.li(3),'Color', [0.8 0.1 0.1], 'LineWidth',1.5, 'LineStyle',':')
%    title('distribution of events across experimental run')
xlim([0,min([DATA.totaltime])/1000])
ylim(get(gca,'YLim'))
h.l4 = line(repmat(linspace(1/6,5/6,5).*min([DATA.totaltime])/1000,2,1),repmat(get(gca,'YLim'),5,1)',...
    'Color', [0.2 0.6 0.2 1],'LineWidth',1.5);
xlabel('time in s')

figs{3} = figure;
set(gcf,'Position',[100 100 600 200],'PaperPositionMode','auto')
subplot(4,1,1)
t.respidx = strcmp([DATA_df.onset_type],'resp');
h.l1 = line([[DATA_df(t.blnkidx).onset_times]./1000; [DATA_df(t.blnkidx).onset_times]./1000],...
    [ones(1,sum(t.blnkidx)); ones(1,sum(t.blnkidx))+0.5], ...
    'Color', [0.2 0.4 1 0.05],'LineWidth',1.5);
hold on
set(gca,'YTickLabels',[],'YTick',[],'XTickLabels',[],'XTick',[],'Box','on')
xlim([0,min([DATA.totaltime])/1000])
ylim([0.95,1.55])
% plot vertical lines
h.l3 = line(repmat(linspace(1/6,5/6,5).*min([DATA.totaltime])/1000,2,1),repmat(get(gca,'YLim'),5,1)',...
    'Color', [0.2 0.6 0.2 1],'LineWidth',1.5);
%     xlabel('time in s')
%     legend([h.l1(1) h.l2(1)], {'finger taps';'blinks'},'Location','southoutside','Orientation','horizontal')
title('distribution of blinks across experimental run')

subplot(4,1,[2 3 4])
h.hist1 = histogram([DATA_df(t.blnkidx).onset_times]./1000,60,'BinLimits',[0,min([DATA.totaltime])/1000],...
    'FaceColor', [0.2 0.4 1]);
hold on
h.li = plot(perm.bincenter,[perm.sh_percentiles_blnk]);
set(h.li(1),'Color', [0.2 0.4 1], 'LineWidth',1.5, 'LineStyle',':')
set(h.li(2),'Color', [0.2 0.4 1], 'LineWidth',1.5, 'LineStyle','-')
set(h.li(3),'Color', [0.2 0.4 1], 'LineWidth',1.5, 'LineStyle',':')
%    title('distribution of events across experimental run')
xlim([0,min([DATA.totaltime])/1000])
ylim(get(gca,'YLim'))
h.l4 = line(repmat(linspace(1/6,5/6,5).*min([DATA.totaltime])/1000,2,1),repmat(get(gca,'YLim'),5,1)',...
    'Color', [0.2 0.6 0.2 1],'LineWidth',1.5);
xlabel('time in s')

%%%%%%%%%%
% histogram of difference between blinks
figs{4} = figure;
set(gcf,'Position',[100 100 600 400],'PaperPositionMode','auto')
subplot(2,1,1)
histogram([DATA.resp_timediffs_intrial]./1000,60,'FaceColor', [0.8 0.1 0.1])
xlabel('time in s')
title('time between consecutive finger taps')
annotation('textbox',[.55 .7 .35 .2],...
    'String',sprintf('Median = %1.3f;\n95%%-interval = [%1.3f, %1.3f]',prctile([DATA.resp_timediffs_intrial]./1000,[50 2.5 97.5])),...
    'EdgeColor','none')

subplot(2,1,2)
histogram([DATA.blink_timediffs_intrial]./1000,60,'FaceColor', [0.2 0.4 1])
xlabel('time in s')
title('time between consecutive blinks')
annotation('textbox',[.55 .21 .35 .2],...
    'String',sprintf('Median = %1.3f;\n95%%-interval = [%1.3f, %1.3f]',prctile([DATA.blink_timediffs_intrial]./1000,[50 2.5 97.5])),...
    'EdgeColor','none')

figs{5} = figure;
set(gcf,'Position',[100 100 600 200],'PaperPositionMode','auto')
histogram([DATA.resp_timediffs_intrial]./1000,60,'FaceColor', [0.8 0.1 0.1])
xlabel('time in s')
title('time between consecutive finger taps')
annotation('textbox',[.55 .6 .35 .2],...
    'String',sprintf('Median = %1.3f;\n95%%-interval = [%1.3f, %1.3f]',prctile([DATA.resp_timediffs_intrial]./1000,[50 2.5 97.5])),...
    'EdgeColor','none')

figs{6} = figure;
set(gcf,'Position',[100 100 600 200],'PaperPositionMode','auto')
histogram([DATA.blink_timediffs_intrial]./1000,60,'FaceColor', [0.2 0.4 1])
xlabel('time in s')
title('time between consecutive blinks')
annotation('textbox',[.55 .6 .35 .2],...
    'String',sprintf('Median = %1.3f;\n95%%-interval = [%1.3f, %1.3f]',prctile([DATA.blink_timediffs_intrial]./1000,[50 2.5 97.5])),...
    'EdgeColor','none')

t.sumdata = [numel(F.Subjects2Use) mean(cellfun('prodofsize', {DATA.resp_samples} )) std(cellfun('prodofsize', {DATA.resp_samples} ))...
     mean(cellfun('prodofsize', {DATA.blink_samples} )) std(cellfun('prodofsize', {DATA.blink_samples} ))];
fprintf(1,['For N = %1.0f subjects\na mean number of %1.3f finger taps (SD = %1.3f) and\n'...
    'a mean number of %1.3f blinks (SD = %1.3f) was recorded.\n'],...
    [t.sumdata])



sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'BlinkResp_time_distribution'; 'Resp_time_distribution'; 'Blink_time_distribution'; ...
    'BlinkResp_timediff_hist'; 'Resp_timediff_hist'; 'Blink_timediff_hist'};
% for i_fig = 1:6
%     print(figs{i_fig}, fullfile(sav.pathout,sav.filenames{i_fig}),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sav.filenames{i_fig}),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sav.filenames{i_fig}),'-depsc2', '-painters','-r300')
% end
% print(fullfile(p.pathout_fft,'figures',sprintf('TOPO_Grand_mean2_%s_%s',pl.dat2plot,p.freqlabel)),'-djpeg','-r300')
% saveas(gcf,fullfile(p.pathout_fft,'figures',sprintf('TOPO_Grand_mean2_%s_%s',pl.dat2plot,p.freqlabel)),'fig')
% print(fullfile(p.pathout_fft,'figures',sprintf('TOPO_Grand_mean2_%s_%s',pl.dat2plot,p.freqlabel)),'-depsc2', '-painters','-r300')

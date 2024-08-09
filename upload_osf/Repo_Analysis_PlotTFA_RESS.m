%% plot previously calculated TFA data

    
clearvars

%% parameters
F.PathIn                = '...osf\eeg_data'; % gabor of 1Hz FWHM [no blinks -3.5 3.5] RESS altered for analysis on Auswerterechner

F.Subjects2Use          = [1 2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20]; % based on trial number for each subject

TFA.baseline            = [-3250 -3000]; % for publication
% TFA.baseline            = [-3250 -2250]; % for revision as asked by a reviewer


pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% loop for subjects
for i_sub=1:numel(F.Subjects2Use)
    %%
    % read in TFA data
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_tfa.mat ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    
    temp.tfa = open(sprintf('%s\\VP%02.0f_exp_tfa.mat',F.PathIn,F.Subjects2Use(i_sub)));
    
    % convert to single
    temp.tfa.TFA.data_induced = single(temp.tfa.TFA.data_induced);
    temp.tfa.TFA.data_evoked = single(temp.tfa.TFA.data_evoked);
    temp.tfa.TFA.data_RESS_induced = single(temp.tfa.TFA.data_RESS_induced);
    temp.tfa.TFA.data_RESS_evoked = single(temp.tfa.TFA.data_RESS_evoked);
    
    % preallocate memory
    if i_sub == 1
        TFA.data_induced = single(nan([size(temp.tfa.TFA.data_induced) numel(F.Subjects2Use)]));
        TFA.data_evoked = TFA.data_induced;
        TFA.data_induced_bc = TFA.data_induced;
        TFA.data_evoked_bc = TFA.data_induced;
        TFA.data_RESS_induced = single(nan([size(temp.tfa.TFA.data_RESS_induced) numel(F.Subjects2Use)]));
        TFA.data_RESS_evoked = TFA.data_RESS_induced ;
        TFA.data_RESS_induced_bc = TFA.data_RESS_induced ;
        TFA.data_RESS_evoked_bc = TFA.data_RESS_induced ;
        TFA.time = temp.tfa.TFA.t;
        TFA.frequency = temp.tfa.TFA.f;
        TFA.electrodes = temp.tfa.TFA.electrodes;
        TFA.electrodes_RESS = temp.tfa.TFA.electrodes_RESS;
        TFA.RESS_map = nan([length(temp.tfa.TFA.RESS.map) numel(F.Subjects2Use)]);
        TFA.RESS_map_signchan = nan([length(temp.tfa.TFA.RESS.map) numel(F.Subjects2Use)]);
        [TFA.RESS_snr_ind, TFA.RESS_snr_evo] = deal(nan(1,numel(F.Subjects2Use)));
        TFA.srate = temp.tfa.TFA.params.srate/2;
    end
    
    % subject by subject
    TFA.data_induced(:,:,:,i_sub) = temp.tfa.TFA.data_induced; % induced data
    TFA.data_induced_bc(:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_induced, ...
        mean(temp.tfa.TFA.data_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
   
    TFA.data_evoked(:,:,:,i_sub) = temp.tfa.TFA.data_evoked; % evoked data
    TFA.data_evoked_bc(:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_evoked, ...
        mean(temp.tfa.TFA.data_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
    
    TFA.data_RESS_induced(:,:,i_sub) = temp.tfa.TFA.data_RESS_induced; % induced data
    TFA.data_RESS_induced_bc(:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_RESS_induced, ...
        mean(temp.tfa.TFA.data_RESS_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
    
    TFA.data_RESS_evoked(:,:,i_sub) = temp.tfa.TFA.data_RESS_evoked; % evoked data
    TFA.data_RESS_evoked_bc(:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_RESS_evoked, ...
        mean(temp.tfa.TFA.data_RESS_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
    
    try TFA.trialnum_raw_induced(i_sub)=numel(temp.tfa.TFA.art.resp_time_n);
    end
    try TFA.trialnum_raw_evoked(i_sub)=numel(temp.tfa.TFA.art.resp_time_sh);
    end
   
    TFA.totaltrials(i_sub)=temp.tfa.TFA.alltrials; % number of total trials
    
    TFA.RESS_map(:,i_sub) = temp.tfa.TFA.RESS.map(:,end);
    TFA.RESS_map_signchan(:,i_sub) = temp.tfa.TFA.RESS.map_signchanged(:,end);
    TFA.RESS_snr_ind(i_sub)=temp.tfa.TFA.RESS.SNR_ind(end);
    TFA.RESS_snr_evo(i_sub)=temp.tfa.TFA.RESS.SNR_evo(end);
    
    clear temp
    
end

%% actual plotting no RESS data
% figure; topoplot([],TFA.electrodes, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];

% plotting parameters
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha
% pl.elec2plot = {'C3';'CP3'};  sav.chan_add = 'MotorSmall';% motor alpha/beta
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.xlims=[-3250 3250]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

pl.flims = [3 35];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);
clear figs fig

pl.fre2index = [8 13.167 15.167 30];

%%%%%%%%%%%%%
% raw induced no RESS
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced(:,:,pl.elec2plot_i,pl.subs2use),3),4));
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));
% t.clims1=[min(min(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))) max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)))];
% t.clims2=[min(min(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))) max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)))];
figs{1}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
subplot(2,1,1)
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for induced activity\n for channel [%s]', vararg2str(pl.elec2plot)), 'FontSize',8)

xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)
cb = colorbar();



subplot(2,1,2)
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for evoked activity\n for channel [%s]', vararg2str(pl.elec2plot)), 'FontSize',8)
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)
cb = colorbar();

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});



% normalized induced no RESS
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[-1 1]*max(max(abs(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))));
figs{2}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
h.sp1=subplot(2,1,1);

imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca,flipud(cbrewer2('RdBu')))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for induced activity; baseline [%1.0f %1.0f]ms\n for channel [%s]',...
    TFA.baseline,vararg2str(pl.elec2plot)),'FontSize',8)
h.cb1=colorbar;
set(h.cb1,'FontSize',8)
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)

t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));


h.sp2=subplot(2,1,2);
% ttest
t.datain = squeeze(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3));
[tt.h tt.p tt.ci tt.stats]=ttest(permute(t.datain, [3 1 2]));
TFA.pvals_induced_bc=squeeze(tt.p);
tt.stats.tstat = squeeze(tt.stats.tstat);



t.data = abs(log10(TFA.pvals_induced_bc(t.find1:t.find2,t.ind1:t.ind2)));
t.clim = [0 max(t.data(:))];
t.pcriterion = abs(log10(0.05));
if max(t.data(:))<t.pcriterion
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time,TFA.frequency,abs(log10(TFA.pvals_induced_bc)),t.clim)
set(gca,'YDir','normal')
title(sprintf('TFA-p-values for induced activity; baseline [%1.0f %1.0f]ms\n for channel [%s]',...
    TFA.baseline,vararg2str(pl.elec2plot)),'FontSize',8)
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
% h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% plot cluster contour lines
if flag_cluster == 1
    hold on
    contour(TFA.time,TFA.frequency,clust.cluster_h_map,[1,1],'Color','c','LineWidth',2)
end

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});



% normalized evoked noRESS
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims2=[-1 1]*max(max(abs(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))));
figs{3}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
h.sp1=subplot(2,1,1);
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
colormap(gca,flipud(cbrewer2('RdBu')))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for evoked activity; baseline [%1.0f %1.0f]ms\n for channel [%s]',...
    TFA.baseline,vararg2str(pl.elec2plot)),'FontSize',8)
% freezeColors

h.cb1=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
set(h.cb1,'FontSize',8)
% h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)

h.sp2=subplot(2,1,2);
for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(nanmean(TFA.data_evoked_bc(i_freq,:,pl.elec2plot_i,pl.subs2use),3))');
    TFA.pvals_induced_bc(i_freq,:)=tt.p;
end
t.data = abs(log10(TFA.pvals_induced_bc(t.find1:t.find2,:)));
t.clim = [0 max(t.data(:))];
t.pcriterion = abs(log10(0.05));
if max(t.data(:))<t.pcriterion
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time,TFA.frequency,abs(log10(TFA.pvals_induced_bc)),t.clim)
set(gca,'YDir','normal')
title(sprintf('TFA-p-values for evoked activity; baseline [%1.0f %1.0f]ms\n for channel [%s]',...
    TFA.baseline,vararg2str(pl.elec2plot)),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
hline(14.16667,'g')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});
 
sav.pathout = '...';
sav.filenames = {'Resp_TFA_Amp_RAW_EvoIndu';'Resp_TFA_Amp_BC_Indu'};
% for i_fig = 1:2
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-depsc2', '-vector','-r300')
% end





%% plot RESS components
pl.xlims=[-3250 3250]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

pl.flims = [3 35];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);
clear figs

% raw induced RESS
pl.data_ind = squeeze(nanmean(TFA.data_RESS_induced(:,:,pl.subs2use),3));
pl.data_evo = squeeze(nanmean(TFA.data_RESS_evoked(:,:,pl.subs2use),3));
t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));
figs{1}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')

subplot(2,1,1)
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca, fake_parula)
set(gca,'YDir','normal')
title(sprintf('tfa RESS component; averaged in frequency domain\nRESS'),'FontSize',8)
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
set(gca,'FontSize',8)

subplot(2,1,2)
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
colormap(gca, fake_parula)
set(gca,'YDir','normal')
title(sprintf('tfa RESS component; averaged in time domain\nRESS'),'FontSize',8)
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
set(gca,'FontSize',8)



% normalized induced RESS
pl.data_ind = squeeze(nanmean(TFA.data_RESS_induced_bc(:,:,pl.subs2use),3));
t.clims1=[-1 1]*max(max(abs(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))));
figs{2}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')

h.sp1=subplot(2,1,1);

% colormap(cmap('R2'))
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca,flipud(cbrewer2('RdBu')))
set(gca,'YDir','normal')
title(sprintf('normalized tfa RESS component;\n averaged in frequency domain; baseline [%1.0f %1.0f]ms', TFA.baseline),'FontSize',8)
h.cb1=colorbar;
set(h.cb1,'FontSize',8)
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
% freezeColors
% h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')
set(gca,'FontSize',8)

h.sp2=subplot(2,1,2);
for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(TFA.data_RESS_induced_bc(i_freq,:,pl.subs2use))');
    TFA.pvals_RESS_induced_bc(i_freq,:)=tt.p;
end
t.data = abs(log10(TFA.pvals_RESS_induced_bc(t.find1:t.find2,:)));
t.clim = [0 max(t.data(:))];
t.pcriterion = abs(log10(0.05));
if max(t.data(:))<t.pcriterion
   t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time,TFA.frequency,t.data,t.clim)
set(gca,'YDir','normal')
title(sprintf('p-values RESS component;\n averaged in frequency domain; baseline [%1.0f %1.0f]ms', TFA.baseline),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
hline(14.16667,'g')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% normalized evoked RESS
pl.data_evo = squeeze(nanmean(TFA.data_RESS_evoked_bc(:,:,pl.subs2use),3));
t.clims2=[-1 1]*max(max(abs(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))));
figs{3}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')

h.sp1=subplot(2,1,1);
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
set(gca,'YDir','normal')
colormap(gca,flipud(cbrewer2('RdBu')))
title(sprintf('normalized tfa RESS component;\n averaged in time domain; baseline [%1.0f %1.0f]ms', TFA.baseline),'FontSize',8)
% freezeColors

h.cb1=colorbar;
set(h.cb1,'FontSize',8)
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)

% h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')
set(gca,'FontSize',8)

h.sp2=subplot(2,1,2);
for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(TFA.data_RESS_evoked_bc(i_freq,:,pl.subs2use))');
    TFA.pvals_RESS_evoked_bc(i_freq,:)=tt.p;
end
t.data = abs(log10(TFA.pvals_RESS_evoked_bc(t.find1:t.find2,:)));
t.clim = [0 max(t.data(:))];
t.pcriterion = abs(log10(0.05));
if max(t.data(:))<t.pcriterion
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time,TFA.frequency,t.data,t.clim)
set(gca,'YDir','normal')
title(sprintf('p-values RESS component;\n averaged in time domain; baseline [%1.0f %1.0f]ms', TFA.baseline),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
% h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% sav.pathout = '...';
% sav.filenames = {'Resp_TFA_Amp_RAW_EvoIndu_RESS';'Resp_TFA_Amp_BC_Indu_RESS';'Resp_TFA_Amp_BC_Evo_RESS'};
% for i_fig = 1:3
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-depsc2', '-painters','-r300')
% end






%% plot all lines of interest into one graphics (separate graphics + single subject data)
clear figs h

% for review
pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 13.16667] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
    [8 13.16667] {'C3';'CP3'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
    [15.16667 30] {'C3';'CP3'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
    };

pl.subs2use = 1:numel(F.Subjects2Use);

pl.xlim = [-3250 3250];

pl.plotnum = [ceil(sqrt(size(pl.parameters,1))) round(sqrt(size(pl.parameters,1)))*3+5];
pl.plpos = reshape(sort(...
    [1:11:prod(pl.plotnum) 2:11:prod(pl.plotnum) 3:11:prod(pl.plotnum) 5:11:prod(pl.plotnum) 6:11:prod(pl.plotnum) 7:11:prod(pl.plotnum)]...
    ),3,[]);
pl.topopos = reshape(sort([8:11:prod(pl.plotnum) 9:11:prod(pl.plotnum) 10:11:prod(pl.plotnum) 11:11:prod(pl.plotnum)]),2,[]);
pl.permut_n = 10000;
% pl.permut_n = 100;


h.pl = [];
figs{1}=figure;
set(gcf,'Position',[100 100 800 550],'PaperPositionMode','auto')
clear pl.tdata2 h.sp h.p h.pm
for i_pl = 1:size(pl.parameters,1)
    h.sp(i_pl,1)=subplot(pl.plotnum(1),pl.plotnum(2),pl.plpos(:,i_pl));
    
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.xlim(2)));
    
    % index electrodes
    switch  pl.parameters{i_pl,7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data=squeeze(nanmean(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.elec2plot_i,pl.subs2use),1),3));',pl.parameters{i_pl,3});
        case 'RESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes_RESS.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.subs2use),1));',pl.parameters{i_pl,3});
        
    end
    eval(com)
    [tt.h tt.p tt.ci tt.stats]=ttest(pl.data');
    % do cluster based permutation
    [cluster_runt_diff, timecourse_runt_diff]=eeg_erpStat_clusterP(pl.data,zeros(size(pl.data)),pl.permut_n,2);
    tt.p(timecourse_runt_diff.h_corr==0)=1;
       
    % actual plot with boundedline
    pl.mdata = mean(pl.data,2)';
    pl.ddata = tt.ci - repmat(mean(pl.data,2)',2,1);
    [h.l(i_pl), h.p(i_pl)] = boundedline(TFA.time(t.ind3:t.ind4),pl.mdata,pl.ddata(1,:));
    set(h.l(i_pl),'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5});
    set(h.p(i_pl),'FaceColor',pl.parameters{i_pl,4},'FaceAlpha',0.25)
    hold on
    
    
    % significant effects?
    pl.tdata2=nan(size(pl.mdata));
    pl.tdata2(tt.p<.05)=1;
    set(gca,'ylim',[-1.1 (1+0.15+0.05*4)]...
        .*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
    pl.maxdata = max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false)));
    pl.y = 1*pl.maxdata*(1+0.05+0.05*2);
    plot(TFA.time(t.ind3:t.ind4),pl.tdata2.*pl.y,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',6)

    % display some effects numerically
    if any(~isnan(pl.tdata2))
        fprintf('results for %s:\n',pl.parameters{i_pl,6})
        % positive cluster
        if ~isempty(cluster_runt_diff.positive.h)
            t.clustidx = find(cluster_runt_diff.positive.h==1);
            for i_clust = 1:numel(t.clustidx)
                % extrct relevant data
                t.clust.time = TFA.time(cluster_runt_diff.positive.samples{t.clustidx(i_clust)}([1 end])+t.ind3-1);
                t.clust.mod_avg = mean(pl.mdata(cluster_runt_diff.positive.samples{t.clustidx(i_clust)}));
                t.clust.mod_sd = std(mean(pl.data(cluster_runt_diff.positive.samples{t.clustidx(i_clust)},:)));
                t.clust.mod_sdwithin = std(pl.mdata(cluster_runt_diff.positive.samples{t.clustidx(i_clust)}));
                t.clust.mod_max = max(pl.mdata(cluster_runt_diff.positive.samples{t.clustidx(i_clust)}));
                t.clust.mod_tmax = max(tt.stats.tstat(cluster_runt_diff.positive.samples{t.clustidx(i_clust)}));
                t.clust.mod_Dmax = t.clust.mod_tmax/sqrt(size(pl.data,2));
                t.clust.mod_tavg = mean(tt.stats.tstat(cluster_runt_diff.positive.samples{t.clustidx(i_clust)}));
                t.clust.mod_Davg = t.clust.mod_tavg/sqrt(size(pl.data,2));
                t.clust.mod_D_cluster = mean(mean(pl.data(cluster_runt_diff.positive.samples{t.clustidx(i_clust)},:))) / ...
                    std(mean(pl.data(cluster_runt_diff.positive.samples{t.clustidx(i_clust)},:))); % as in publication
                t.clust.MonteCarloP = cluster_runt_diff.positive.p(t.clustidx(t.clustidx(i_clust)));

%                 % check cluster
%                 [t.clust.tt.h t.clust.tt.p t.clust.tt.ci t.clust.tt.stats]= ...
%                     ttest(mean(pl.data(cluster_runt_diff.positive.samples{t.clustidx(i_clust)},:)));
%                 t.clust.tt.stats.tstat/sqrt(size(pl.data,2));

                % fprintf
                fprintf(['positive cluster %1.0f: from %1.0f to %1.0f ms; avg mod = %1.3f; sd mod = %1.3f; ' ...
                    'max mod = %1.3f; max t = %1.3f; max D = %1.3f; cluster p = %1.3f; cluster D = %1.3f\n'], ...
                    i_clust, t.clust.time, t.clust.mod_avg, t.clust.mod_sd, t.clust.mod_max, t.clust.mod_tmax, t.clust.mod_Dmax, t.clust.MonteCarloP, t.clust.mod_D_cluster)
            end
        end

        % negative cluster
        if ~isempty(cluster_runt_diff.negative.h)
            t.clustidx = find(cluster_runt_diff.negative.h==1);
            for i_clust = 1:numel(t.clustidx)
                % extrct relevant data
                t.clust.time = TFA.time(cluster_runt_diff.negative.samples{t.clustidx(i_clust)}([1 end])+t.ind3-1);
                t.clust.mod_avg = mean(pl.mdata(cluster_runt_diff.negative.samples{t.clustidx(i_clust)}));
                t.clust.mod_sd = std(mean(pl.data(cluster_runt_diff.negative.samples{t.clustidx(i_clust)},:)));
                t.clust.mod_sdwithin = std(pl.mdata(cluster_runt_diff.negative.samples{t.clustidx(i_clust)}));
                t.clust.mod_max = min(pl.mdata(cluster_runt_diff.negative.samples{t.clustidx(i_clust)}));
                t.clust.mod_tmax = min(tt.stats.tstat(cluster_runt_diff.negative.samples{t.clustidx(i_clust)}));
                t.clust.mod_Dmax = t.clust.mod_tmax/sqrt(size(pl.data,2));
                t.clust.mod_tavg = mean(tt.stats.tstat(cluster_runt_diff.negative.samples{t.clustidx(i_clust)}));
                t.clust.mod_Davg = t.clust.mod_tavg/sqrt(size(pl.data,2));
                t.clust.mod_D_cluster = mean(mean(pl.data(cluster_runt_diff.negative.samples{t.clustidx(i_clust)},:))) / ...
                    std(mean(pl.data(cluster_runt_diff.negative.samples{t.clustidx(i_clust)},:))); % as in publication
                t.clust.MonteCarloP = cluster_runt_diff.negative.p(t.clustidx(t.clustidx(i_clust)));

                % fprintf
                fprintf(['negative cluster %1.0f: from %1.0f to %1.0f ms; avg mod = %1.3f; sd mod = %1.3f; ' ...
                    'max mod = %1.3f; max t = %1.3f; max D = %1.3f; cluster p = %1.3f; cluster D = %1.3f\n'], ...
                    i_clust, t.clust.time, t.clust.mod_avg, t.clust.mod_sd, t.clust.mod_max, t.clust.mod_tmax, t.clust.mod_Dmax, t.clust.MonteCarloP, t.clust.mod_D_cluster)
            end
        end
    end
    
    
    % graphics
    xlim(pl.xlim)
    title(pl.parameters{i_pl,6},'FontSize', 8)
    grid on
    box on
    hline(0,'k')
    vline(TFA.baseline,'k')
    xlabel('time in ms')
    ylabel('amplitude modulation in %')
    
end
% %%%%%%% vers1 %start%
% set(gca,'ylim',[-1.1 1.1].*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
% %%%%%%% vers1 %end%

for i_pl = 1:size(pl.parameters,1)
    subplot(pl.plotnum(1),pl.plotnum(2),pl.topopos(:,i_pl))
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
%     topoplot([],TFA.electrodes(1:64),'style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',2});
    topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
        'emarker2',{find(pl.elec2plot_i),'o',pl.parameters{i_pl,4},3,1});
    title(sprintf('%s\n[%1.1f %1.1f]Hz',pl.parameters{i_pl,6},pl.parameters{i_pl,1}),'Color',pl.parameters{i_pl,4})
end

sav.pathout = '...';
sav.filenames = {'Resp_AllSignals_Amp_Timecourse_sep_v2_9_rev01'};
% for i_fig = 1:1
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-depsc2', '-vector','-r300')
% end



%% extract data for time window to do some explorative analyses
clear figs h



% for review
pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 13.16667] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
%     [8 13.16667] {'C3';'CP3'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
%     [15.16667 30] {'C3';'CP3'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
    };

pl.subs2use = 1:numel(F.Subjects2Use);

pl.time2use = [-1014 229]; % visual alpha cluster


clearvars cluster_runt_diff timecourse_runt_diff tt 
pl.data = []; 
for i_pl = 1:size(pl.parameters,1)
    
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.time2use(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.time2use(2)));
    
        
    % index electrodes
    switch  pl.parameters{i_pl,7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data{i_pl}=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.elec2plot_i,pl.subs2use),[1,2,3]));',pl.parameters{i_pl,3});
            
        case 'RESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes_RESS.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data{i_pl}=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.subs2use),[1,2]));',pl.parameters{i_pl,3});
            
    end
    eval(com)
    [tt(i_pl).h tt(i_pl).p tt(i_pl).ci tt(i_pl).stats]=ttest(pl.data{i_pl}');
    % do cluster based permutation
    
end

% test
[correl.rho,correl.pval] = corr(pl.data{1},pl.data{2}, 'type','Spearman');

% save data
dataout_t = table(F.Subjects2Use(pl.subs2use)',pl.data{1},pl.data{2},'VariableNames',{'particpant','SSVEP','vis.alpha'});


t.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_volmov\data_in';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile
% writetable(dataout_t,fullfile(t.path,sprintf('GaborModulations_%1.0fto%1.0fms%s.csv',pl.time2use,t.datestr)),'Delimiter',';')

%% correlational analysis add on
[correl.rho correl.pval] = corr(pl.data{1},pl.data{3}(:,:)', 'type','Spearman');
fdr(correl.pval)

 figure;
 topoplot(correl.rho, TFA.electrodes, ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',[-1 1], ...
                'colormap',flipud(cbrewer2('RdBu')),'whitebk','on');
        title(sprintf('significant electrodes in cluster'))

    colorbar












    


%% topoplot
% pl.time2plot=[500 2000];
% pl.time2plot=[-1000 0];
pl.time2plot=[-500 300];
% pl.time2plot=[-2000 500];
% pl.time2plot=[-2000 750];
[t.t t.ind1]=min(abs(TFA.time-pl.time2plot(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.time2plot(2)));

% pl.freq2plot=[14.16667 14.16667];
% pl.freq2plot=[10 14];
% pl.freq2plot=[15 25];
% pl.freq2plot=[20 30];
pl.freq2plot=[8 12];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.data_raw=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_raw(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_bc=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_bc(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));

pl.subs2use = 1:numel(F.Subjects2Use);

pl.titleadd = {'induced';'evoked'};

figure;
% loop for evoked/induced
for i_dat = 1:2
    h.s1=subplot(2,3,1+3*(i_dat-1));
    pl.data=squeeze(nanmean(pl.data_raw(:,:,i_dat),2));
    topoplot( pl.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(pl.data)],'conv','on','colormap',fake_parula);
    title(sprintf('raw %s amp %1.2f to %1.2f Hz;\n%1.0f to %1.0f ms',pl.titleadd{i_dat},pl.freq2plot, pl.time2plot))
%     freezeColors
    h.c1 = colorbar;
%     h.c1=cbfreeze(h.c1);
    
    h.s2=subplot(2,3,2+3*(i_dat-1));
    pl.data=squeeze(nanmean(pl.data_bc(:,:,i_dat),2));
    topoplot( pl.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits','absmax','conv','on','colormap',colormap(gca,flipud(cbrewer2('RdBu'))));
    title(sprintf('normalized %s amp %1.2f to %1.2f Hz;\n%1.0f to %1.0f ms',pl.titleadd{i_dat},pl.freq2plot, pl.time2plot))
%     freezeColors
    h.c2 = colorbar;
%     h.c2=cbfreeze(h.c2);
    
    [tt.h tt.p tt.ci tt.stats]=ttest(pl.data_bc(:,:,i_dat)');
    t.data = abs(log10(tt.p));
    t.clim = [0 max(t.data(:))];
    t.pcriterion = abs(log10(0.05));
    if max(t.data(:))<t.pcriterion
        % temp.colormap = repmat([0.5 0.5 0.5],100,1);
        t.colormap = repmat(linspace(1,0.3,1000)',1,3);
    else
        t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
        % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
        t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
    end
    h.s2=subplot(2,3,3+3*(i_dat-1));
    topoplot( t.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',t.clim,'conv','on','colormap',t.colormap);
    title(sprintf('p-values %s %1.2f to %1.2f Hz;\n%1.0f to %1.0f ms',pl.titleadd{i_dat},pl.freq2plot, pl.time2plot))
%     freezeColors
    h.cb2 = colorbar;
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))
%     h.cb2=cbfreeze(h.cb2);
    
end

%% topoplot of only single topos
clear figs

% specific topos for effects [revision
% pl.time2plot=[-1014 229]; pl.freq2plot=[8 13.167]; %vis. alpha
% pl.time2plot=[-1412 494]; pl.freq2plot=[8 13.167]; %motor. alpha
pl.time2plot=[1314 2705]; pl.freq2plot=[8 13.167]; %motor. alpha
pl.time2plot=[-1333.98 228.52]; pl.freq2plot=[15.167 30]; %motor. beta
pl.time2plot=[400.39 3250]; pl.freq2plot=[15.167 30]; %motor. beta



[t.t t.ind1]=min(abs(TFA.time-pl.time2plot(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.time2plot(2)));

[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.elec2plot = {};
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.subs2use = 1:numel(F.Subjects2Use);            

pl.data_raw=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_raw(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_bc=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_bc(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));

% colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))



figs{1}=figure;
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
pl.data=squeeze(nanmean(pl.data_raw(:,:,1),2));
topoplot( pl.data, TFA.electrodes(1:64), ...
    'shading', 'flat', 'numcontour', 0, 'maplimits',[0 max(pl.data)], 'colormap', fake_parula,...
    'conv','on','whitebk','on','emarker2', {find(pl.elec2plot_i),'o','m',8});

title(sprintf('raw; %1.2f to %1.2f Hz; %1.0f to %1.0f ms',pl.freq2plot, pl.time2plot))
%     freezeColors
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])



figs{2}=figure;
pl.data=squeeze(nanmean(pl.data_bc(:,:,1),2));
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
topoplot( pl.data, TFA.electrodes(1:64),...
    'shading', 'flat', 'numcontour', 0,'maplimits','absmax', 'colormap', flipud(cbrewer2('RdBu')),...
    'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
title(sprintf('bc; %1.2f to %1.2f Hz; %1.0f to %1.0f ms',pl.freq2plot, pl.time2plot))
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])


% sav.pathout = '...';
sav.filenames = {'TOPO_raw_effects';'TOPO_bc_effects'};
% for i_fig = 1:2
%     print(figs{i_fig}, fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-dpng','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-depsc2', '-vector','-r300')
% end


%% plot RESS.topography

pl.data = mean(TFA.RESS_map_signchan(:,:),2);
figs{1}=figure;
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
topoplot( pl.data, TFA.electrodes(1:64),'shading', 'flat', 'numcontour', 0,'maplimits','absmax',...
    'whitebk','on','conv','on','colormap',fake_parula);
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])

sav.pathout = '...';
sav.filenames = {'TOPO_RESS'};
% for i_fig = 1:1
% %     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-depsc2', '-painters','-r300')
% end

%% plot electrode head
% pl.elec2plot = {'C3';'CP3'};
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};


figure;
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
% topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on');
topoplot([],TFA.electrodes(1:64),'whitebk','on','style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',8});



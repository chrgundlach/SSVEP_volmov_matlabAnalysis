%% plot previously calculated TFA data

    
clearvars
%% parameters
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS';
% F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks';
% F.PathIn                = 'F:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks2'; % latest used
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD'; % gabor of 1Hz FWHM
% F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD2'; % gabor of 0.5Hz FWHM
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD3'; % gabor of 1Hz FWHM [no blinks -3 3]
F.PathIn                = 'C:\Matlab_user\Christopher\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD4'; % gabor of 1Hz FWHM [no blinks -3.5 3.5]
% F.Subjects2Use          = [1:4];
F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % based on trial number for each subject
% F.Subjects2Use          = [1 2 3 4];

% TFA.baseline            = [-3000 -2750];
% TFA.baseline            = [-2750 -2500];
% TFA.baseline            = [-2250 -2000];
TFA.baseline            = [-3250 -3000];
% TFA.baseline            = [-4000 -3500];
% TFA.baseline            = [-4500 -4000];



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
%     TFA.data_induced_bc(:,:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_induced, ...
%         mean(temp.tfa.TFA.data_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
    TFA.data_induced_bc(:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_induced, ...
        mean(temp.tfa.TFA.data_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
   
    TFA.data_evoked(:,:,:,i_sub) = temp.tfa.TFA.data_evoked; % evoked data
%     TFA.data_evoked_bc(:,:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_evoked, ...
%         mean(temp.tfa.TFA.data_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
    TFA.data_evoked_bc(:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_evoked, ...
        mean(temp.tfa.TFA.data_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
    
    TFA.data_RESS_induced(:,:,i_sub) = temp.tfa.TFA.data_RESS_induced; % induced data
%     TFA.data_RESS_induced_bc(:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_RESS_induced, ...
%         mean(temp.tfa.TFA.data_RESS_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
    TFA.data_RESS_induced_bc(:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_RESS_induced, ...
        mean(temp.tfa.TFA.data_RESS_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
    
    TFA.data_RESS_evoked(:,:,i_sub) = temp.tfa.TFA.data_RESS_evoked; % evoked data
%     TFA.data_RESS_evoked_bc(:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_RESS_evoked, ...
%         mean(temp.tfa.TFA.data_RESS_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
    TFA.data_RESS_evoked_bc(:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_RESS_evoked, ...
        mean(temp.tfa.TFA.data_RESS_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
    
    try TFA.trialnum_raw_induced(i_sub)=numel(temp.tfa.TFA.art.resp_time_n);
    end
    try TFA.trialnum_raw_evoked(i_sub)=numel(temp.tfa.TFA.art.resp_time_sh);
    end
    try TFA.trialnum_discarded_induced(i_sub)=numel(temp.tfa.TFA.art.SCADS_Trials2Del_n);
    end
    try TFA.trialnum_discarded_evoked(i_sub)=numel(temp.tfa.TFA.art.SCADS_Trials2Del_sh);
    end
    TFA.totaltrials(i_sub)=temp.tfa.TFA.alltrials; % number of total trials
    %TFA.trials_analyzed_evoked(:,i_sub)=(cellfun(@(x) numel(x),temp.tfa.TFA.trials_evoked)./temp.tfa.TFA.alltrials)*100; % percentage analysed trials
    
    TFA.RESS_map(:,i_sub) = temp.tfa.TFA.RESS.map(:,end);
    TFA.RESS_map_signchan(:,i_sub) = temp.tfa.TFA.RESS.map_signchanged(:,end);
    TFA.RESS_snr_ind(i_sub)=temp.tfa.TFA.RESS.SNR_ind(end);
    TFA.RESS_snr_evo(i_sub)=temp.tfa.TFA.RESS.SNR_evo(end);
    
    clear temp
    
end

%% actual plotting no RESS data
% figure; topoplot([],TFA.electrodes, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);


% plotting parameters
% pl.elec2plot = {'Oz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; % vis alpha
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
% pl.elec2plot = {'POz'}; % vis alpha II
% pl.elec2plot = {'C3';'CP3'}; % motor alpha/beta
pl.elec2plot = {'C3';'CP3';'C5';'CP5'};  sav.chan_add = 'MotorLarge';% motor alpha/beta II
% pl.elec2plot = {'F1';'F2';'Fz'}; sav.chan_add = 'Frontal1'; % frontal
% pl.elec2plot = {'FP1';'FP2'}; sav.chan_add = 'Frontal2'; % frontal
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; sav.chan_add = 'Parietal'; % parietal
% pl.elec2plot = {'F5';'F3';'F1';'FC3'}; sav.chan_add = 'PreMotor'; % pre-motor
% pl.elec2plot = {'FC3';'FC3';'FC1';'FCz';'Fz';'F5';'F3';'F1';'AF7';'AF3';'AFz';'FP1';'FPz'};  sav.chan_add = 'CentroFrontalLeft'; % pre-motor
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% pl.xlims=[-2500 2500]; % index time 2 plot
pl.xlims=[-3250 3250]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

% pl.flims = TFA.frequency([1 end]); % index frequency 2 plot
pl.flims = [3 TFA.frequency(end)]; % index frequency 2 plot
% pl.flims = [4 TFA.frequency(end)]; % index frequency 2 plot
% pl.flims = [4 35];
% pl.flims = [3 35];
% pl.flims = [8 44];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
clear figs fig


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
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

% if logscale

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
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();
% cb.Ruler.Scale = 'log';
% cb.Ruler.MinorTick = 'on';

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
pl.colmap = [[linspace(0,0.9,128)' linspace(0,0.9,128)' linspace(1,0.9,128)']; [linspace(0.9,1,128)' linspace(0.9,0,128)' linspace(0.9,0,128)']];
% pl.colmap = [[zeros(64,1) linspace(1,0,64)' ones(64,1)];...
%     [linspace(0,1,64)' linspace(0,1,64)' linspace(1,1,64)'];...
%     [linspace(1,1,64)' linspace(1,0,64)' linspace(1,0,64)'];...
%     [ones(64,1) linspace(0,1,64)' zeros(64,1)]];
% colormap(cmap('R2'))
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
% colormap(gca,pl.colmap); % cbrewer('div','RdBu',256,'l')
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
% freezeColors
% h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')
set(gca,'FontSize',8)

t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));


h.sp2=subplot(2,1,2);
for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(nanmean(TFA.data_induced_bc(i_freq,:,pl.elec2plot_i,pl.subs2use),3))');
    TFA.pvals_induced_bc(i_freq,:)=tt.p;
end
t.data = abs(log10(TFA.pvals_induced_bc(t.find1:t.find2,t.ind1:t.ind2)));
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
imagesc(TFA.time,TFA.frequency,abs(log10(TFA.pvals_induced_bc)),t.clim)
set(gca,'YDir','normal')
title(sprintf('TFA-p-values for induced activity; baseline [%1.0f %1.0f]ms\n for channel [%s]',...
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
% h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

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
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
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
% h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'Resp_TFA_Amp_RAW_EvoIndu';'Resp_TFA_Amp_BC_Indu';'Resp_TFA_Amp_BC_Evo'};
for i_fig = 1:3
    print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'fig')
    print(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-depsc2', '-painters','-r300')
end



%% plot only single fields

% plotting parameters
% pl.elec2plot = {'Oz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; % vis alpha
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; % vis alpha II
% pl.elec2plot = {'POz'}; % vis alpha II
% pl.elec2plot = {'POz'}; % vis alpha IV
% pl.elec2plot = {'C3';'CP3'}; % motor alpha/beta
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; % motor alpha/beta II
% pl.elec2plot = {'F1';'F2';'Fz'}; % frontal
% pl.elec2plot = {'Fp1';'Fp2'}; % frontal
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% pl.xlims=[-2500 2500]; % index time 2 plot
pl.xlims=[-3500 3500]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

pl.flims = TFA.frequency([1 end]); % index frequency 2 plot
pl.flims = [4 35];
% pl.flims = [8 44];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2



%%%%%%%%%%%%%
% raw induced no RESS
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced(:,:,pl.elec2plot_i,pl.subs2use),3),4));
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));
figure;
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca, fake_parula)
set(gca,'YDir','normal')
title(sprintf('tfa noRESS, for channel %s;\n averaged in frequency domain', vararg2str(pl.elec2plot)), 'FontSize',8)
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
set(gca,'FontSize',8)


% raw normalized no RESS
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[-1 1]*max(max(abs(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))));
figure;
colormap(gca,flipud(cbrewer2('RdBu')))
% colormap(gca,cmap('R1'))
% colormap(gca,cmap('D9'))
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
set(gca,'YDir','normal')
title(sprintf('normalized tfa noRESS, for channel %s;\n averaged in frequency domain; baseline [%1.0f %1.0f]ms',...
    vararg2str(pl.elec2plot),TFA.baseline),'FontSize',8)
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

for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(nanmean(TFA.data_induced_bc(i_freq,:,pl.elec2plot_i,pl.subs2use),3))');
    TFA.pvals_induced_bc(i_freq,:)=tt.p;
end
pl.alpha = ones(size(TFA.pvals_induced_bc)).*0.85;
pl.alpha(TFA.pvals_induced_bc<.05) = 1;
% alpha(pl.alpha)


%% statistically testing differences in TFA images between signals from different sources
% figure; topoplot([],TFA.electrodes, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);

% plotting parameters
pl.parameters = [];
% vis alpha II large
pl.parameters(1).elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};
pl.parameters(1).name = 'VisualLarge';
% motor alpha/beta II
pl.parameters(end+1).elec2plot = {'C3';'CP3';'C5';'CP5'};
pl.parameters(end).name = 'MotorLarge';
% % parietal
% pl.parameters(end+1).elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'};
% pl.parameters(end).name = 'Parietal';

for i_pl = 1:numel(pl.parameters)
    pl.parameters(i_pl).elec2plot_i = logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.parameters(i_pl).elec2plot, 'UniformOutput',false)),1));
end

% pl.xlims=[-2500 2500]; % index time 2 plot
pl.xlims=[-3250 3250]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

% pl.flims = TFA.frequency([1 end]); % index frequency 2 plot
% pl.flims = [3 TFA.frequency(end)]; % index frequency 2 plot
pl.flims = [4 TFA.frequency(end)]; % index frequency 2 plot
% pl.flims = [4 35];
% pl.flims = [3 35];
% pl.flims = [8 44];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.baseline = TFA.baseline; %[-3250 -3000];
pl.baseline_i = dsearchn(TFA.time', pl.baseline');

pl.subs2use = 1:numel(F.Subjects2Use);
clear figs fig


%extract data first
pl.testdata = nan(numel(t.find1:t.find2),numel(t.ind1:t.ind2),numel(pl.parameters),numel(F.Subjects2Use),2);
% freq, time, electrode cluster, subjects, [baseline time]
for i_pl = 1:numel(pl.parameters)
    pl.testdata(:,:,i_pl,:,1)= mean(TFA.data_induced(t.find1:t.find2,t.ind1:t.ind2,pl.parameters(i_pl).elec2plot_i,pl.subs2use),3);
    pl.testdata(:,:,i_pl,:,2)= repmat(...
       mean(mean(TFA.data_induced(t.find1:t.find2,pl.baseline_i(1):pl.baseline_i(2),pl.parameters(i_pl).elec2plot_i,pl.subs2use),2),3),...
       [1,size(pl.testdata,2),1,1]);
end

% do liner mixed modelling for each frequency/time combination

prog = 0; i_up = 0;TFA.lme=[];
fprintf('doing %1.0f lme_fits; progress: %5.2f%%\n',size(pl.testdata,1)*size(pl.testdata,2),prog)
for i_freq = 1:size(pl.testdata,1)
    for i_time = 1:size(pl.testdata,2)
        prog = ( 100*(i_up/(size(pl.testdata,1)*size(pl.testdata,2))));
        fprintf(1,'\b\b\b\b\b\b%5.2f%%',prog);
        % data in?
        t.data = squeeze(pl.testdata(i_freq,i_time,:,:,:));
        t.fct1 = repmat({pl.parameters(:).name}',[1 size(t.data,2) size(t.data,3)]);
        t.fct2 = repmat(pl.subs2use,[size(t.data,1) 1 size(t.data,3)]);
        t.fct3 = permute(repmat([1;2], [1 size(t.data,1) size(t.data,2)]),[2 3 1]);
        
        t.tbldata = table(t.data(:),t.fct1(:),t.fct2(:),t.fct3(:),...
            'VariableNames',{'amplitude','source','participant','time'});
        % any systematic change of amplitude due to source and time?
        
        % fit lme
        t.lme = fitlme(t.tbldata,'amplitude~source*time+(1|participant)');
        TFA.lme.F(i_freq,i_time).source = t.lme.anova.FStat(2);
        TFA.lme.p(i_freq,i_time).source = t.lme.anova.pValue(2);
        TFA.lme.F(i_freq,i_time).time = t.lme.anova.FStat(3);
        TFA.lme.p(i_freq,i_time).time = t.lme.anova.pValue(3);
        TFA.lme.F(i_freq,i_time).sourceXtime = t.lme.anova.FStat(4);
        TFA.lme.p(i_freq,i_time).sourceXtime = t.lme.anova.pValue(4);
        %         % post-hoc coefficient tests
        %         [t.pVal,t.F,t.DF1,t.DF2] = coefTest(t.lme,[0 0 0 1 0 0])
        i_up = i_up+1;
    end
end

%%
%%%%%% plot results
pl.pldata = mean(diff(pl.testdata,1,5),4);
pl.pldata = mean((((pl.testdata(:,:,:,:,1)./pl.testdata(:,:,:,:,2))-1).*100),4);
pl.clim = [-1 1].*repmat(max(abs(pl.pldata(:))),size(pl.testdata,3),1);
pl.clim = [-1 1].*squeeze(max(max(abs(pl.pldata))));
figure;
set(gcf,'Position',[100 100 1600 800],'PaperPositionMode','auto')
for i_pl = 1:size(pl.testdata,3)
    
    subplot(3,3,i_pl)
    imagesc(TFA.time(t.ind1:t.ind2),TFA.frequency(t.find1:t.find2),pl.pldata(:,:,i_pl),pl.clim(i_pl,:))
%     colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    colorbar
    title(sprintf('TFA-amplitude difference to [%1.0f %1.0f]ms\n for channel [%s]',pl.baseline, vararg2str(pl.parameters(i_pl).elec2plot)), 'FontSize',8)
    xlabel('time in ms')
    ylabel('difference to baseline')
    ylim(TFA.frequency([t.find1 t.find2]))
    xlim(TFA.time([t.ind1 t.ind2]))
    hold on; plot(TFA.time([t.ind1 t.ind2]),[14.16667 14.1677],'r')
    
    % draw topograpgy
    h.a1 = axes('position',[(0.04 + (0.85/3)*(i_pl-1)) 0.77 0.08 0.08],'Visible','off');
    topoplot(pl.parameters(i_pl).elec2plot_i ,TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
        'emarker2',{find(pl.parameters(i_pl).elec2plot_i),'o','r',4,1});
end


% draw F-values
subplot(3,3,4)
pl.tdata = reshape([TFA.lme.F(:,:).source],size(TFA.lme.F));
imagesc(TFA.time(t.ind1:t.ind2),TFA.frequency(t.find1:t.find2),pl.tdata,[0 max(pl.tdata(:))])
%     colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
colorbar
title(sprintf('F-values SOURCE'), 'FontSize',8)
xlabel('time in ms')
ylabel('difference to baseline')
ylim(TFA.frequency([t.find1 t.find2]))
xlim(TFA.time([t.ind1 t.ind2]))
hold on; plot(TFA.time([t.ind1 t.ind2]),[14.16667 14.1677],'r')

subplot(3,3,5)
pl.tdata = reshape([TFA.lme.F(:,:).time],size(TFA.lme.F));
imagesc(TFA.time(t.ind1:t.ind2),TFA.frequency(t.find1:t.find2),pl.tdata,[0 max(pl.tdata(:))])
%     colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
colorbar
title(sprintf('F-values TIME'), 'FontSize',8)
xlabel('time in ms')
ylabel('difference to baseline')
ylim(TFA.frequency([t.find1 t.find2]))
xlim(TFA.time([t.ind1 t.ind2]))
hold on; plot(TFA.time([t.ind1 t.ind2]),[14.16667 14.1677],'r')

subplot(3,3,6)
pl.tdata = reshape([TFA.lme.F(:,:).sourceXtime],size(TFA.lme.F));
imagesc(TFA.time(t.ind1:t.ind2),TFA.frequency(t.find1:t.find2),pl.tdata,[0 max(pl.tdata(:))])
%     colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
colorbar
title(sprintf('F-values SOURCE X TIME'), 'FontSize',8)
xlabel('time in ms')
ylabel('difference to baseline')
ylim(TFA.frequency([t.find1 t.find2]))
xlim(TFA.time([t.ind1 t.ind2]))
hold on; plot(TFA.time([t.ind1 t.ind2]),[14.16667 14.1677],'r')


% draw p-values
subplot(3,3,7)
pl.tdata = abs(log10(reshape([TFA.lme.p(:,:).source],size(TFA.lme.p))));
t.clim = [0 max(pl.tdata(:))];
t.pcriterion = abs(log10(0.05));
if max(pl.tdata(:))<t.pcriterion
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(pl.tdata(:)))*1000);
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time(t.ind1:t.ind2),TFA.frequency(t.find1:t.find2),pl.tdata,t.clim)
set(gca,'YDir','normal')
title(sprintf('p-values'),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
% h.cb2=cbfreeze(h.cb2);
colormap(gca,t.colormap)


subplot(3,3,8)
pl.tdata = abs(log10(reshape([TFA.lme.p(:,:).time],size(TFA.lme.p))));
t.clim = [0 max(pl.tdata(:))];
t.pcriterion = abs(log10(0.05));
if max(pl.tdata(:))<t.pcriterion
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(pl.tdata(:)))*1000);
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time(t.ind1:t.ind2),TFA.frequency(t.find1:t.find2),pl.tdata,t.clim)
set(gca,'YDir','normal')
title(sprintf('p-values'),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
% h.cb2=cbfreeze(h.cb2);
colormap(gca,t.colormap)

subplot(3,3,9)
pl.tdata = abs(log10(reshape([TFA.lme.p(:,:).sourceXtime],size(TFA.lme.p))));
t.clim = [0 max(pl.tdata(:))];
t.pcriterion = abs(log10(0.05));
if max(pl.tdata(:))<t.pcriterion
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(pl.tdata(:)))*1000);
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time(t.ind1:t.ind2),TFA.frequency(t.find1:t.find2),pl.tdata,t.clim)
set(gca,'YDir','normal')
title(sprintf('p-values'),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
% h.cb2=cbfreeze(h.cb2);
colormap(gca,t.colormap)


%% plot RESS components
pl.xlims=[-3250 3250]; % index time 2 plot
% pl.xlims=[-3500 3500]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

% pl.flims = TFA.frequency([1 end]); % index frequency 2 plot
% pl.flims = [4 44];
% pl.flims = [4 35];
pl.flims = [3 35];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
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
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
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
% h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% normalized evoked RESS
pl.data_evo = squeeze(nanmean(TFA.data_RESS_evoked_bc(:,:,pl.subs2use),3));
t.clims2=[-1 1]*max(max(abs(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))));
figs{3}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')

h.sp1=subplot(2,1,1);
% colormap(cmap('R2'))
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

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'Resp_TFA_Amp_RAW_EvoIndu_RESS';'Resp_TFA_Amp_BC_Indu_RESS';'Resp_TFA_Amp_BC_Evo_RESS'};
for i_fig = 1:3
    print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'fig')
    print(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-depsc2', '-painters','-r300')
end

%% plot timecourse of frequency
% pl.freq2plot=[14.16667 14.16667];
% pl.freq2plot=[10 14];
% pl.freq2plot=[15 25];
% pl.freq2plot=[20 30];
pl.freq2plot=[8 14];
% pl.freq2plot=[15 30];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(end)));

% pl.elec2plot = {'Oz';'POz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; % vis alpha I
% pl.elec2plot = {'POz'}; % vis alpha II
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; % vis alpha III
% pl.elec2plot = {'POz';'Pz';'P1';'P2'}; % vis alpha II
% pl.elec2plot = {'C3';'CP3'}; % motor alpha/beta
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; % motor alpha/beta II
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.xlims=[-2750 3000]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

pl.subs2use = F.Subjects2Use;

pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));

% ttest for each point
[tt.h pl.pdata_ind tt.ci tt.stats]=ttest(squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,:,pl.elec2plot_i,:),1),3))');
[tt.h pl.pdata_evo tt.ci tt.stats]=ttest(squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,:,pl.elec2plot_i,:),1),3))');

figure;
subplot(2,1,1)
plot(TFA.time,mean(pl.data_ind(t.ind3:t.ind4,:),1),'k')
hold on
try
   pl.tdata=mean(pl.data_ind(t.ind3:t.ind4,:),1);
   pl.tdata(pl.pdata_ind>=.05)=nan;
   plot(TFA.time,pl.tdata,'b','LineWidth',2)
   try
       pl.tdata(pl.pdata_ind>=.01)=nan;
       plot(TFA.time,pl.tdata,'c','LineWidth',2)
       try
           pl.tdata(pl.pdata_ind>=.001)=nan;
           plot(TFA.time,pl.tdata,'r','LineWidth',2)
           legend({'amplitude';'p<.05';'p<.01';'p<.001'})
       catch
           legend({'amplitude';'p<.05';'p<.01'})
       end
   catch
       legend({'amplitude';'p<.05'})
   end
catch 
    legend('amplitude')
end
ylim([-1.2 1.2]*max(abs(mean(pl.data_ind(t.ind3:t.ind4,t.ind1:t.ind2),1))))
xlabel('time in ms')
ylabel('amplitude in \muV')
xlim(pl.xlims)
title(sprintf('induced timecourse signal %1.2f to %1.2f Hz\n channel %s; baseline [%1.0f %1.0f]ms',...
    pl.freq2plot, vararg2str(pl.elec2plot),TFA.baseline))
grid on
hline(0,'k')

subplot(2,1,2)
plot(TFA.time,mean(pl.data_evo(t.ind3:t.ind4,:),1),'k')
hold on
try
   pl.tdata=mean(pl.data_evo(t.ind3:t.ind4,:),1);
   pl.tdata(pl.pdata_evo>=.05)=nan;
   plot(TFA.time,pl.tdata,'b','LineWidth',2)
   try
       pl.tdata(pl.pdata_evo>=.01)=nan;
       plot(TFA.time,pl.tdata,'c','LineWidth',2)
       try
           pl.tdata(pl.pdata_evo>=.001)=nan;
           plot(TFA.time,pl.tdata,'r','LineWidth',2)
           legend({'amplitude';'p<.05';'p<.01';'p<.001'})
       catch
           legend({'amplitude';'p<.05';'p<.01'})
       end
   catch
       legend({'amplitude';'p<.05'})
   end
catch 
    legend('amplitude')
end
ylim([-1.2 1.2]*max(abs(mean(pl.data_evo(t.ind3:t.ind4,t.ind1:t.ind2),1))))
xlabel('time in ms')
ylabel('amplitude in \muV')
xlim(pl.xlims)
title(sprintf('evoked timecourse signal %1.2f to %1.2f Hz\n channel %s; baseline [%1.0f %1.0f]ms',...
    pl.freq2plot, vararg2str(pl.elec2plot),TFA.baseline))
grid on
hline(0,'k')

%% plot all lines of interest into one graphics

% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_induced_bc' [0.6 0.6 0.6] 2 sprintf('vis ind\nSSVEP') 'RESS';...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
%     [8 12] {'PO4';'O2';'PO8';'P8';'P10';'I2'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
%     [18 30] {'PO4';'O2';'PO8';'P8';'P10';'I2'} 'data_induced_bc' [1 0 1] 1 'vis beta' 'noRESS';...
%     [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
%     };
pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_induced_bc' [0.6 0.6 0.6] 2 sprintf('vis ind\nSSVEP') 'RESS';...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
    [18 30] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 1] 1 'vis beta' 'noRESS';...
    [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
    };
% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_induced_bc' [0.6 0.6 0.6] 2 sprintf('vis ind\nSSVEP') 'RESS';...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
%     [8 12] {'POz'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
%     [18 30] {'POz'} 'data_induced_bc' [1 0 1] 1 'vis beta' 'noRESS';...
%     [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
%     };

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
% 
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % subjects discarded due to low trial number version 2


% pl.xlim = [-3250 3250];
pl.xlim = [-4000 4000];
pl.xlim = [-3500 3500];

pl.toponum = [ceil(sqrt(size(pl.parameters,1))) round(sqrt(size(pl.parameters,1)))];
pl.plpos = repmat([1:pl.toponum(2)*3],pl.toponum(1),1)+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2)*3);
pl.topopos = repmat([pl.toponum(2)*3+1:pl.toponum(2)*4],pl.toponum(1),1)'+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2))';

h.pl = [];
figure;
subplot(pl.toponum(1),pl.toponum(2)*4,pl.plpos(:))
for i_pl = 1:size(pl.parameters,1)
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
       
    % actual plot
    h.pl(i_pl)=plot(TFA.time(t.ind3:t.ind4),mean(pl.data,2),'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5});
    hold on
    
    % add significant portions
    pl.tdata = nan(size(pl.data));
    pl.t2data = mean(pl.data,2);
    try
        pl.tdata(tt.p<.05)=pl.t2data(tt.p<.05);
    end
    plot(TFA.time(t.ind3:t.ind4),pl.tdata,'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5}+2)
end
xlim(pl.xlim)
set(gca,'ylim',[-1.1 1.1].*max(cellfun(@(x) max(abs(x)), get(get(gca,'Children'),'YData'))))
legend(h.pl, pl.parameters(:,6),'location','NorthWest','FontSize', 8)
grid on
hline(0,'k')
xlabel('time in ms')
ylabel('amplitude')

for i_pl = 1:size(pl.parameters,1)
    subplot(pl.toponum(1),pl.toponum(2)*4,pl.topopos(i_pl))
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
    topoplot([],TFA.electrodes(1:64),'style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',2});
    title(pl.parameters{i_pl,6},'Color',pl.parameters{i_pl,4})
end

%% plot all lines of interest into one graphics (bounded line)
clear figs h
% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_induced_bc' [0.6 0.6 0.6] 2 sprintf('vis ind\nSSVEP') 'RESS';...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
%     [8 12] {'PO4';'O2';'PO8';'P8';'P10';'I2'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
%     [18 30] {'PO4';'O2';'PO8';'P8';'P10';'I2'} 'data_induced_bc' [1 0 1] 1 'vis beta' 'noRESS';...
%     [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
%     };
pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_induced_bc' [0.6 0.6 0.6] 2 sprintf('vis ind\nSSVEP') 'RESS';...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [18 30] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 1] 2 'vis beta' 'noRESS';...
    [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
    };
% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_induced_bc' [0.6 0.6 0.6] 2 sprintf('vis ind\nSSVEP') 'RESS';...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
%     [8 12] {'POz'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
%     [18 30] {'POz'} 'data_induced_bc' [1 0 1] 1 'vis beta' 'noRESS';...
%     [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
%     };
% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
%     [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
%     [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
%     };
% 
% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
%     [8 12] {'PO4';'O2';'PO8';'P8';'P10';'I2'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
%     [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
%     };
% 
pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
    };

pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
    [16 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0.8 0 0.8] 2 'parietal beta' 'noRESS';...
    };

pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0.8 0 0.8] 2 'parietal beta' 'noRESS';...
    };

% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
%     [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
%     [15 21] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [0.4 0.4 0.4] 2 'vis beta' 'noRESS';...
%     [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
%     [16 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0.8 0 0.8] 2 'parietal beta' 'noRESS';...
%     [8 12] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0 0.8 0.8] 2 'parietal alpha' 'noRESS';...
%     };




pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
% 
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % subjects discarded due to low trial number version 2

% pl.xlim = [-3250 3250];
pl.xlim = [-4000 4000];
pl.xlim = [-3500 3500];

pl.toponum = [ceil(sqrt(size(pl.parameters,1))) round(sqrt(size(pl.parameters,1)))];
pl.plpos = repmat([1:pl.toponum(2)*3],pl.toponum(1),1)+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2)*3);
pl.topopos = repmat([pl.toponum(2)*3+1:pl.toponum(2)*4],pl.toponum(1),1)'+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2))';
pl.permut_n = 10000;
% pl.permut_n = 100;


h.pl = [];
figs{1}=figure;
set(gcf,'Position',[100 100 800 350],'PaperPositionMode','auto')
h.sp(1)=subplot(pl.toponum(1),pl.toponum(2)*4,pl.plpos(:));
clear pl.tdata2
for i_pl = 1:size(pl.parameters,1)
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
%     [cluster_runt_diff, timecourse_runt_diff]=eeg_erpStat_clusterP(pl.data,zeros(size(pl.data)),pl.permut_n,2);
%     tt.p(timecourse_runt_diff.h_corr==0)=1;
       
    % actual plot with boundedline
    pl.mdata = mean(pl.data,2)';
    pl.ddata = tt.ci - repmat(mean(pl.data,2)',2,1);
    %pl.ddata = tt.ci;
    [h.l(i_pl), h.p(i_pl)] = boundedline(TFA.time(t.ind3:t.ind4),pl.mdata,pl.ddata(1,:));
    set(h.l(i_pl),'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5});
    set(h.p(i_pl),'FaceColor',pl.parameters{i_pl,4},'FaceAlpha',0.25)
    hold on
    
%     %%%%%%% vers1 %start%
%     % add significant portions
%     pl.tdata = nan(size(pl.mdata));
%     pl.t2data = mean(pl.data,2);
%     try
%         pl.tdata(tt.p<.05)=pl.t2data(tt.p<.05);
%     end
%     plot(TFA.time(t.ind3:t.ind4),pl.tdata,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',pl.parameters{i_pl,5}+3)
%     %%%%%%% vers1 %end%
    
    pl.tdata2(i_pl,:)=nan(size(pl.mdata));
    pl.tdata2(i_pl,tt.p<.05)=1; 
end
xlim(pl.xlim)
% %%%%%%% vers1 %start%
% set(gca,'ylim',[-1.1 1.1].*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
% %%%%%%% vers1 %end%

%%%%%%% vers2 %start%
% plot lines for significant effects below
set(gca,'ylim',[-(1+0.15+0.05*size(pl.parameters,1)) 1.1].*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
pl.maxdata = max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false)));
for i_pl = 1:size(pl.parameters,1)
    pl.y = -1*pl.maxdata*(1+0.05+0.05*i_pl);
    
%     plot(TFA.time(t.ind3:t.ind4),repmat(pl.y,1,numel(TFA.time(t.ind3:t.ind4))),'Color','k','LineWidth',0.1)
    plot(TFA.time(t.ind3:t.ind4),pl.tdata2(i_pl,:).*pl.y,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',2)
    
end
%%%%%%% vers2 %end%

legend(h.l, pl.parameters(:,6),'location','NorthWest','FontSize', 8)
grid on
box on
hline(0,'k')
xlabel('time in ms')
ylabel('amplitude modulation in %')

for i_pl = 1:size(pl.parameters,1)
    subplot(pl.toponum(1),pl.toponum(2)*4,pl.topopos(i_pl))
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
%     topoplot([],TFA.electrodes(1:64),'style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',2});
    topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
        'emarker2',{find(pl.elec2plot_i),'o','r',3,1});
    title(sprintf('%s\n[%1.1f %1.1f]Hz',pl.parameters{i_pl,6},pl.parameters{i_pl,1}),'Color',pl.parameters{i_pl,4})
end

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_5'};
sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2_clustcorr_5'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2_clustcorr_5b'};
for i_fig = 1:1
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-depsc2', '-painters','-r300')
end

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
% pl.time2plot=[500 2500];
% pl.time2plot=[500 2000];
% pl.time2plot=[-1000 0];
% pl.time2plot=[-3000 -2750];
% pl.time2plot=[-500 300];
% pl.time2plot=[-2000 500];
% pl.time2plot=[-1500 0];
pl.time2plot=[1500 2500];
% pl.time2plot=[-3500 3500];
[t.t t.ind1]=min(abs(TFA.time-pl.time2plot(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.time2plot(2)));

% pl.freq2plot=[14.16667 14.16667];
pl.freq2plot=[10 14];
% pl.freq2plot=[16 23];
% pl.freq2plot=[15 23];
% pl.freq2plot=[18 30];
pl.freq2plot=[8 12];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.elec2plot = {'C3';'CP3'};
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};
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
% topoplot( pl.data, TFA.electrodes(1:64), ...
%     'shading', 'flat', 'numcontour', 0, 'maplimits',[0 max(pl.data)], 'colormap', flipud(cbrewer2('greys')),...
%     'conv','on','whitebk','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
% topoplot( pl.data, TFA.electrodes(1:64), ...
%     'shading', 'flat', 'numcontour', 0, 'maplimits','minmax', 'colormap', cbrewer2('greys'),...
%     'conv','on','whitebk','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
title(sprintf('raw; %1.2f to %1.2f Hz; %1.0f to %1.0f ms',pl.freq2plot, pl.time2plot))
%     freezeColors
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])
% set(h.cb,'Position',[t.pos(1)+0.1 t.pos(2) t.pos(3) t.pos(4)])
%     h.c1=cbfreeze(h.c1);


figs{2}=figure;
pl.data=squeeze(nanmean(pl.data_bc(:,:,1),2));
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
topoplot( pl.data, TFA.electrodes(1:64),...
    'shading', 'flat', 'numcontour', 0,'maplimits','absmax', 'colormap', flipud(cbrewer2('RdBu')),...
    'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
% topoplot( pl.data, TFA.electrodes(1:64),...
%     'shading', 'flat', 'numcontour', 0,'maplimits',[-10 10], 'colormap', flipud(cbrewer2('RdBu')),...
%     'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
title(sprintf('bc; %1.2f to %1.2f Hz; %1.0f to %1.0f ms',pl.freq2plot, pl.time2plot))
%     freezeColors
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])


sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'TOPO_raw';'TOPO_bc'};
for i_fig = 1:2
%     print(figs{i_fig}, fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-dpng','-r300')
    print(figs{i_fig}, fullfile(sav.pathout,...
        sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,...
        sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'fig')
    print(figs{i_fig},fullfile(sav.pathout,...
        sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-depsc2', '-painters','-r300')
end

%% topoplot differences...rudimentary
clear figs
% time windows of to be subtracted data
pl.time2plot = {[500 2000];[500 2000]};
% pl.time2plot = {[-1000 0];[-1000 0]};
pl.time2plot_i=cell2mat(cellfun(@(x) dsearchn(TFA.time',x'),pl.time2plot,'UniformOutput',false)');

pl.freq2plot={[16 23];[18 30]};
% pl.freq2plot={[8 12];[10 14]};
pl.freq2plot_i=cell2mat(cellfun(@(x) dsearchn(TFA.frequency',x'),pl.freq2plot,'UniformOutput',false)');

pl.elec2plot = {'C3';'CP3'};
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};
pl.elec2plot = {};
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.subs2use = 1:numel(F.Subjects2Use);            

pl.data_raw=squeeze(nanmean(nanmean(TFA.data_induced(pl.freq2plot_i(1,1):pl.freq2plot_i(2,1),pl.time2plot_i(1,1):pl.time2plot_i(2,1),:,pl.subs2use),1),2))-...
    squeeze(nanmean(nanmean(TFA.data_induced(pl.freq2plot_i(1,2):pl.freq2plot_i(2,2),pl.time2plot_i(1,2):pl.time2plot_i(2,2),:,pl.subs2use),1),2));
pl.data_raw(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(pl.freq2plot_i(1,1):pl.freq2plot_i(2,1),pl.time2plot_i(1,1):pl.time2plot_i(2,1),:,pl.subs2use),1),2))-...
    squeeze(nanmean(nanmean(TFA.data_evoked(pl.freq2plot_i(1,2):pl.freq2plot_i(2,2),pl.time2plot_i(1,2):pl.time2plot_i(2,2),:,pl.subs2use),1),2));
pl.data_bc=squeeze(nanmean(nanmean(TFA.data_induced_bc(pl.freq2plot_i(1,1):pl.freq2plot_i(2,1),pl.time2plot_i(1,1):pl.time2plot_i(2,1),:,pl.subs2use),1),2))-...
    squeeze(nanmean(nanmean(TFA.data_induced_bc(pl.freq2plot_i(1,2):pl.freq2plot_i(2,2),pl.time2plot_i(1,2):pl.time2plot_i(2,2),:,pl.subs2use),1),2));
pl.data_bc(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(pl.freq2plot_i(1,1):pl.freq2plot_i(2,1),pl.time2plot_i(1,1):pl.time2plot_i(2,1),:,pl.subs2use),1),2))-...
    squeeze(nanmean(nanmean(TFA.data_evoked_bc(pl.freq2plot_i(1,2):pl.freq2plot_i(2,2),pl.time2plot_i(1,2):pl.time2plot_i(2,2),:,pl.subs2use),1),2));
% colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))



figs{1}=figure;
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
pl.data=squeeze(nanmean(pl.data_raw(:,:,1),2));
h.tp=topoplot( pl.data, TFA.electrodes(1:64),...
    'shading', 'flat', 'numcontour', 0,'maplimits','absmax', 'colormap', flipud(cbrewer2('RdBu')),...
    'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
t.pos = get(gca,'Position');
set(gca,'Position',[t.pos(1) t.pos(2)-0.1 t.pos(3) t.pos(4) ])
title(sprintf('raw diff; [%1.1f:%1.1f]-[%1.1f:%1.1f] Hz\n [%1.0f:%1.0f]-[%1.0f:%1.0f] ms',...
    pl.freq2plot{1},pl.freq2plot{2}, pl.time2plot{1},pl.time2plot{1}))
%     freezeColors
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])
% set(h.cb,'Position',[t.pos(1)+0.1 t.pos(2) t.pos(3) t.pos(4)])
%     h.c1=cbfreeze(h.c1);


figs{2}=figure;
pl.data=squeeze(nanmean(pl.data_bc(:,:,1),2));
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
h.tp=topoplot( pl.data, TFA.electrodes(1:64),...
    'shading', 'flat', 'numcontour', 0,'maplimits','absmax', 'colormap', flipud(cbrewer2('RdBu')),...
    'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
h.tp=topoplot( pl.data, TFA.electrodes(1:64),...
    'shading', 'flat', 'numcontour', 0,'maplimits',[-5 5], 'colormap', flipud(cbrewer2('RdBu')),...
    'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
t.pos = get(gca,'Position');
set(gca,'Position',[t.pos(1) t.pos(2)-0.1 t.pos(3) t.pos(4) ])
title(sprintf('bc diff; [%1.1f:%1.1f]-[%1.1f:%1.1f] Hz\n [%1.0f:%1.0f]-[%1.0f:%1.0f] ms',...
    pl.freq2plot{1},pl.freq2plot{2}, pl.time2plot{1},pl.time2plot{1}))
%     freezeColors
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])


sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'TOPO_raw';'TOPO_bc'};
% for i_fig = 1:2
% %     print(figs{i_fig}, fullfile(sav.pathout,...
% %         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-depsc2', '-painters','-r300')
% end

%% plot RESS.topography

% pl.data = mean(TFA.RESS_map(:,:),2);
% pl.data = mean(TFA.RESS_map(:,:).*-1,2);
pl.data = mean(TFA.RESS_map_signchan(:,:),2);
figs{1}=figure;
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
topoplot( pl.data, TFA.electrodes(1:64),'shading', 'flat', 'numcontour', 0,'maplimits','absmax',...
    'whitebk','on','conv','on','colormap',fake_parula);
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'TOPO_RESS'};
for i_fig = 1:1
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-dpng','-r300')
    print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'fig')
    print(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-depsc2', '-painters','-r300')
end

%% plot electrode head
% pl.elec2plot = {'C3';'CP3'};
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};


figure;
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on');
% topoplot([],TFA.electrodes(1:64),'whitebk','on','style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',8});

%% plot topographies across time


% define parameters
t.p_time = [100 100 -3000 3000]; % width step min max
t.p_time = [100 100 -3000 200]; % width step min max
t.p_time = [100 100 -3000 500]; % width step min max
% t.p_time = [100 100 -3000 3200]; % width step min max
t.p_time = [100 100 -2500 3000]; % width step min max
t.posScale = 1.1;

% pl.freq2plot=[14.16667 14.16667];
pl.freq2plot=[10 14];
% pl.freq2plot=[16 23];
% pl.freq2plot=[15 19];
% pl.freq2plot=[8 12];
% pl.freq2plot=[15 30];
% pl.freq2plot=[18 30];
% pl.freq2plot=[2 4];
% pl.freq2plot=[3 8];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.data=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,3)=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,4)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
t.conlabel={'raw induced';'raw evoked';'normalized induced';'normalized evoked'};
t.con_lims = [1 1 0 0]; % 1 = 0 to max; % 0 = minmax

t.time=[];
t.timedot=[];
for i_st = 1:floor((t.p_time(4)-t.p_time(1)-t.p_time(3))/t.p_time(2))+1
    t.time(i_st,:)=t.p_time(3)+t.p_time(2)*(i_st-1)+[0 t.p_time(1)];
    [t.t t.timedot(i_st,1)]=min(abs(TFA.time-t.time(end,1)));
    [t.t t.timedot(i_st,2)]=min(abs(TFA.time-t.time(end,2)));
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.data(t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),1))';
end


for i_fig = 1:size(plotdata,1)
    h.fig(i_fig)=figure;
    if t.con_lims(i_fig) == 1
        t.lims = [0 1]*max(max(abs(plotdata(i_fig,:,:))));
    else
        t.lims = [-1 1]*max(max(abs(plotdata(i_fig,:,:))));
    end
    for i_spl = 1:size(plotdata,3)
        h.sp(i_spl)=subplot(t.row,t.col,i_spl);
        if i_fig < 3
            topoplot( plotdata(i_fig,:,i_spl), TFA.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims,'electrodes','off','colormap',fake_parula);
        else
            topoplot( plotdata(i_fig,:,i_spl), TFA.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims,'electrodes','off','colormap',flipud(cbrewer2('RdBu')));
        end
        title(sprintf('[%1.0f %1.0f]',t.time(i_spl,1),t.time(i_spl,2)),'FontSize',8)
        t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((t.posScale-1)/2)) t.pos(3:4).*t.posScale])
    end
    h.sp(i_spl+1)=subplot(t.row,t.col,i_spl+2);
    topoplot( [], TFA.electrodes(1:64),  ...
        'style','blank');
    title(sprintf('%s\n[%1.2f %1.2f]Hz',t.conlabel{i_fig},pl.freq2plot),'FontSize',8)
    t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((t.posScale-1)/2)) t.pos(3:4).*t.posScale])
    
    t.pos2 = get(h.sp(i_spl+1),'Position');
    t.pos3 = get(h.sp(i_spl+1),'OuterPosition');
    h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
    if i_fig < 3
        colormap(gca, fake_parula)
    else
        colormap(gca,flipud(cbrewer2('RdBu')))
    end
    caxis(t.lims);
    h.c3 = colorbar();
    t.pos4 = get(h.c3,'Position');
    set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3)/2 t.pos4(4)*(2/3)])
end
axcopy(h.fig(i_fig))

%% plot topographies across time [movie]


% define parameters
t.p_time = [100 -3000 3000]; % width step min max

% pl.freq2plot=[14.16667 14.16667];
pl.freq2plot=[10 14];
% pl.freq2plot=[16 25];
% pl.freq2plot=[8 12];
% pl.freq2plot=[15 30];
% pl.freq2plot=[18 30];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.data=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,3)=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,4)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
t.conlabel={'raw induced';'raw evoked';'normalized induced';'normalized evoked'};
t.con_lims = [1 1 0 0]; % 1 = 0 to max; % 0 = minmax


t.timedot=[];
t.t_lims = dsearchn(TFA.time', t.p_time(2:3)');
t.t_steps = t.t_lims(1):t.t_lims(2);
for i_st = 1:numel(t.t_steps)
    % index time windows to be averaged
    t.timedot(i_st,:)=dsearchn(TFA.time',TFA.time(t.t_steps(i_st))+[-1; 1]*t.p_time(1)/2);
end

% create data
plotdata=nan(numel(t.t_lims(1):t.t_lims(2)),size(pl.data,2),size(pl.data,3));
for i_pl = 1:size(t.timedot,1)
    plotdata(i_pl,:,:)=squeeze(mean(pl.data(t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),1));
end

% create movie
for i_dat = 3
    figure;
    if t.con_lims(i_dat) == 1
        pl.clims = [0 1]*max(max(abs(plotdata(:,:,i_dat))));
    else
        pl.clims = [-1 1]*max(max(abs(plotdata(:,:,i_dat))));
    end
    if i_dat < 3
        [mov.Movie,mov.Colormap] = eegmovie(plotdata(:,:,i_dat)',TFA.srate,TFA.electrodes,...
            'topoplotopt',{'maplimits',pl.clims,'numcontour', 0, 'conv','on', 'shading', 'flat','colormap',fake_parula},...
            'time','on','startsec',-3);
    else
        [mov.Movie,mov.Colormap] = eegmovie(plotdata(:,:,i_dat)',TFA.srate,TFA.electrodes,...
            'topoplotopt',{'maplimits',pl.clims,'numcontour', 0, 'conv','on', 'shading', 'flat','colormap',flipud(cbrewer2('RdBu'))},...
            'time','on','startsec',-3);
    end
end
figure; seemovie(mov.Movie,-5,mov.Colormap);


%% plot topographies across time [save single images]
pl.dat2plot = 3; %1=induced raw; 2=evoked raw; 3=induced bc; 4=evoked bc

% sav.pathout = 'C:\Users\psy05cvd\Pictures\gifski\SSVEP_volmov\upper_alpha_topo';
% sav.pathout = 'C:\Users\psy05cvd\Pictures\gifski\SSVEP_volmov\lower_alpha_topo';
sav.pathout = 'C:\Users\psy05cvd\Pictures\gifski\SSVEP_volmov\lower_beta_topo';
% sav.pathout = 'C:\Users\psy05cvd\Pictures\gifski\SSVEP_volmov\higher_beta_topo';


% define parameters
t.p_time = [100 100 -3000 3000]; % width step min max
t.p_time = [100 100 -3000 200]; % width step min max
t.p_time = [100 100 -3000 500]; % width step min max
% t.p_time = [100 100 -3000 3200]; % width step min max
t.p_time = [100 50 -3000 3000]; % width step min max
t.posScale = 1.1;

% pl.freq2plot=[14.16667 14.16667];
% pl.freq2plot=[10 14];
% pl.freq2plot=[16 23];
pl.freq2plot=[15 23];
% pl.freq2plot=[15 19];
% pl.freq2plot=[8 12];
% pl.freq2plot=[15 30];
% pl.freq2plot=[18 30];
% pl.freq2plot=[2 4];
% pl.freq2plot=[3 8];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.data=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,3)=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,4)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
t.conlabel={'raw_induced';'raw_evoked';'normalized_induced';'normalized_evoked'};
t.con_lims = [1 1 0 0]; % 1 = 0 to max; % 0 = minmax

t.time=[];
t.timedot=[];
for i_st = 1:floor((t.p_time(4)-t.p_time(1)-t.p_time(3))/t.p_time(2))+1
    t.time(i_st,:)=t.p_time(3)+t.p_time(2)*(i_st-1)+[0 t.p_time(1)];
    [t.t t.timedot(i_st,1)]=min(abs(TFA.time-t.time(end,1)));
    [t.t t.timedot(i_st,2)]=min(abs(TFA.time-t.time(end,2)));
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.data(t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),1))';
end


fig1=figure;
set(gcf,'Position',[100 100 200 230],'PaperPositionMode','auto')

if t.con_lims(pl.dat2plot) == 1
    t.lims = [0 1]*max(max(abs(plotdata(pl.dat2plot,:,:))));
else
    t.lims = [-1 1]*max(max(abs(plotdata(pl.dat2plot,:,:))));
    t.lims = [-1 1]*10;
end
for i_spl = 1:size(plotdata,3)
    clf(fig1,'reset')
    h.sp1=subplot(4,1,1:3);
    if pl.dat2plot < 3
        topoplot( plotdata(pl.dat2plot,:,i_spl), TFA.electrodes(1:64), ...
            'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims,'electrodes','off',...
            'colormap',fake_parula,'whitebk','on');
    else
        topoplot( plotdata(pl.dat2plot,:,i_spl), TFA.electrodes(1:64), ...
            'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims,'electrodes','off',...
            'colormap',flipud(cbrewer2('RdBu')),'whitebk','on');
    end
    h.cb = colorbar;
    
    %
    h.sp2=subplot(4,1,4);
    plot([min(t.time(:)) max(t.time(:))], [0 0],'Color','k','LineWidth',1)
    hold on
    line([0;0],[-1; 1],'Color',[0 0 0],'LineWidth', 1)
    line([-3000 -2000 -1000 1000 2000 3000;-3000 -2000 -1000 1000 2000 3000],...
        [-0.5 -0.5 -0.5 -0.5 -0.5 -0.5; 0.5  0.5  0.5  0.5  0.5  0.5],'Color',[0 0 0],'LineWidth', 1)
    line([t.time(i_spl,1) t.time(i_spl,2) t.time(i_spl,1) t.time(i_spl,1); t.time(i_spl,2) t.time(i_spl,2) t.time(i_spl,2) t.time(i_spl,1)],...
        [-1 -1 1 -1; -1 1 1 1],'Color',[0.8 0 0],'LineWidth', 2)
    xlim([min(t.time(:)) max(t.time(:))])
    ylim([-1.5 1.5])
    %     set(gca,'YColor','none','XAxisLocation','Origin','Box','Off')
    set(gca,'YColor','none','XColor',[1 1 1],'Box','Off','XTick',[])
    xlabel(sprintf('[%1.0f %1.0f]ms',t.time(i_spl,:)),'Color','k')
    title(sprintf('[%1.0f %1.0f]ms | [%1.0f %1.0f]Hz',[min(t.time(:)) max(t.time(:))], pl.freq2plot))
    
    t.pos = get(h.sp1,'Position');
%     set(h.sp1,'Position',[t.pos(1:2)-(t.pos(3:4).*((t.posScale-1)/2)) t.pos(3:4).*t.posScale])
%     set(h.sp1,'Position',[t.pos(1:2)-(t.pos(3:4).*((t.posScale-1)*2)) t.pos(3:4).*t.posScale])
    set(h.sp1,'Position',[t.pos(1)-(t.pos(3).*((t.posScale-1)*2.5)) t.pos(2)-(t.pos(4).*((t.posScale-1)/2)) t.pos(3:4).*t.posScale])
    t.pos2 = get(h.cb,'Position');
    set(h.cb,'Position',[t.pos2(1) t.pos2(2)+(1/6)*t.pos2(4) t.pos2(3) t.pos2(4)*2/3 ])
    
    % save file
    print(fig1, ...
        fullfile(sav.pathout,sprintf('%s_%1.0f_%1.0fHz_%04.0f',t.conlabel{pl.dat2plot},pl.freq2plot,i_spl)),...
        '-dpng','-r300')
end


%% calculate correlations between amplitudes for frequencies and timewindows of interest
pl.parameters = {...
    [8 12] [-1000 0] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'vis alpha'...
    [8 12] [-1000 0] {'C3';'CP3'} 'motor alpha';...
    [8 12] [-1000 0] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'vis alpha'...
    [10 14] [-1000 0] {'C3';'CP3'} 'motor alpha';...
    [10 14] [-1000 0] {'C3';'CP3'} 'motor alpha'...
    [18 30] [-1000 0] {'C3';'CP3'} 'motor beta';...
    [10 14] [-1000 0] {'C3';'CP3'} 'motor alpha'...
    [10 14] [1500 2500] {'C3';'CP3'} 'motor alpha';...
    [18 30] [-1000 0] {'C3';'CP3'} 'motor beta'...
    [18 30] [1500 2500] {'C3';'CP3'} 'motor beta';...
    };

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
% 
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % subjects discarded due to low trial number version 2


h.pl = [];
for i_pl = 1:size(pl.parameters,1)
    figure;
    % index frequencies
    t.ind1 = dsearchn(TFA.frequency',pl.parameters{i_pl,1}');
    t.ind2 = dsearchn(TFA.frequency',pl.parameters{i_pl,5}');
   
    % index time
    t.ind3 = dsearchn(TFA.time',pl.parameters{i_pl,2}');
    t.ind4 = dsearchn(TFA.time',pl.parameters{i_pl,6}');
    
    % index electrodes
    pl.elec2plot_i1=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,3}, 'UniformOutput',false)),1));
    pl.elec2plot_i2=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,7}, 'UniformOutput',false)),1));
    
    % get xdata
    pl.xdata = squeeze(mean(mean(mean(TFA.data_induced_bc(t.ind1(1):t.ind1(2),t.ind3(1):t.ind3(2),pl.elec2plot_i1,:),1),2),3));
    
    % get ydata
    pl.ydata = squeeze(mean(mean(mean(TFA.data_induced_bc(t.ind2(1):t.ind2(2),t.ind4(1):t.ind4(2),pl.elec2plot_i2,:),1),2),3));
    
    % plotting
    
    h.sc1=scatter(pl.xdata , pl.ydata ,'ko');
    hold on
    
    % xlims
    xlim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
    ylim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
    % correlation
    [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
    text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
        sprintf('R=%1.3f p=%1.4f\n',r.R(2),r.P(2)), ...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
    % add fitted line
    r.fit = polyfit(pl.xdata, pl.ydata,1);
    t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
    
    t.pos=get(gca,'Position');
    
    % title and labels
    title(sprintf(...
        '%s | %1.1f to %1.1f Hz | %1.0f to %1.0f ms | %s X\n%s | %1.1f to %1.1f Hz | %1.0f to %1.0f ms | %s',...
        pl.parameters{i_pl,4}, pl.parameters{i_pl,1}, pl.parameters{i_pl,2}, vararg2str(pl.parameters{i_pl,3}),...
        pl.parameters{i_pl,8}, pl.parameters{i_pl,5}, pl.parameters{i_pl,6}, vararg2str(pl.parameters{i_pl,7})...
        ),'FontSize',8)
    xlabel(pl.parameters{i_pl,4},'FontSize',8)
    ylabel(pl.parameters{i_pl,8},'FontSize',8)
    
    grid on
    box on
    
    
end
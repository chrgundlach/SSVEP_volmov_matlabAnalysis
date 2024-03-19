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

TFA.baseline            = [-3500 -3000];


pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% loop for subjects
for i_sub=1:numel(F.Subjects2Use)
    %%
    % read in TFA data
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_EOG.mat ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    
    temp.eog = open(sprintf('%s\\VP%02.0f_exp_EOG.mat',F.PathIn,F.Subjects2Use(i_sub)));
    
    % preallocate memory
    if i_sub == 1
        EOG.data.ext_blink_n_mean   = nan(temp.eog.EOG.EEG_ext_blink_ep_n.nbchan,temp.eog.EOG.EEG_ext_blink_ep_n.pnts,numel(F.Subjects2Use));
        EOG.data.ext_blink_sh_mean  = EOG.data.ext_blink_n_mean;
        EOG.data.ext_blink_n_std    = EOG.data.ext_blink_n_mean;
        EOG.data.ext_blink_sh_std   = EOG.data.ext_blink_n_mean;
        EOG.time.ext_blink          = temp.eog.EOG.EEG_ext_blink_ep_n.times;
        EOG.electrodes.ext_blink    = temp.eog.EOG.EEG_ext_blink_ep_n.chanlocs;
        
        EOG.data.ext_resp_n_mean    = nan(temp.eog.EOG.EEG_ext_resp_ep_n.nbchan, temp.eog.EOG.EEG_ext_resp_ep_n.pnts, numel(F.Subjects2Use));
        EOG.data.ext_resp_sh_mean   = EOG.data.ext_resp_n_mean;
        EOG.data.ext_resp_n_std     = EOG.data.ext_resp_n_mean;
        EOG.data.ext_resp_sh_std    = EOG.data.ext_resp_n_mean;
        EOG.time.ext_resp           = temp.eog.EOG.EEG_ext_resp_ep_n.times;
        EOG.electrodes.ext_resp     = temp.eog.EOG.EEG_ext_resp_ep_n.chanlocs;
        
        EOG.data.eeg_blink_n_mean   = nan(temp.eog.EOG.EEG_exp_blink_ep_n_ress.nbchan, temp.eog.EOG.EEG_exp_blink_ep_n_ress.pnts, numel(F.Subjects2Use));
        EOG.data.eeg_blink_sh_mean  = EOG.data.eeg_blink_n_mean;
        EOG.data.eeg_blink_n_std    = EOG.data.eeg_blink_n_mean;
        EOG.data.eeg_blink_sh_std   = EOG.data.eeg_blink_n_mean;
        EOG.time.eeg_blink          = temp.eog.EOG.EEG_exp_blink_ep_n_ress.times;
        EOG.electrodes.eeg_blink    = temp.eog.EOG.EEG_exp_blink_ep_n_ress.chanlocs;
        
        EOG.RESS_map = nan([length(temp.eog.EOG.RESS.map) numel(F.Subjects2Use)]);
        [EOG.RESS_snr_ind, EOG.RESS_snr_evo] = deal(nan(1,numel(F.Subjects2Use)));
        
        EOG.blinkmat_n = [];
        EOG.blinkmat_sh = [];
        
        EOG.srate = temp.eog.EOG.EEG_ext_resp_ep_n.srate;
        
        EOG.surr_resp = temp.eog.EOG.EEG_ext_resp_ep_sh;
        EOG.surr_blink = temp.eog.EOG.EEG_exp_blink_ep_sh_ress;
        
        % TFA
        EOG.data.TFA_blink_noRESS_ind = nan([size(temp.eog.EOG.TFA.data_induced),numel(F.Subjects2Use)]);
        EOG.data.TFA_blink_noRESS_evo = nan([size(temp.eog.EOG.TFA.data_induced),numel(F.Subjects2Use)]);
        EOG.data.TFA_blink_noRESS_ind_bc = nan([size(temp.eog.EOG.TFA.data_induced),numel(F.Subjects2Use)]);
        EOG.data.TFA_blink_noRESS_evo_bc = nan([size(temp.eog.EOG.TFA.data_induced),numel(F.Subjects2Use)]);
        EOG.data.TFA_blink_RESS_ind = nan([1 size(temp.eog.EOG.TFA.data_RESS_evoked),numel(F.Subjects2Use)]);
        EOG.data.TFA_blink_RESS_evo = nan([1 size(temp.eog.EOG.TFA.data_RESS_evoked),numel(F.Subjects2Use)]);
        EOG.data.TFA_blink_RESS_ind_bc = nan([1 size(temp.eog.EOG.TFA.data_RESS_evoked),numel(F.Subjects2Use)]);
        EOG.data.TFA_blink_RESS_evo_bc = nan([1 size(temp.eog.EOG.TFA.data_RESS_evoked),numel(F.Subjects2Use)]);
        EOG.data.TFA_electrodes = temp.eog.EOG.TFA.electrodes;
        EOG.data.TFA_params = temp.eog.EOG.TFA.params;
        EOG.data.TFA_times = temp.eog.EOG.TFA.t;
        EOG.data.TFA_freqs = temp.eog.EOG.TFA.f;
    end
    
    % extract data
    EOG.data.ext_blink_n_mean(:,:,i_sub)    = mean(temp.eog.EOG.EEG_ext_blink_ep_n.data,3);
    EOG.data.ext_blink_sh_mean(:,:,i_sub)   = mean(temp.eog.EOG.EEG_ext_blink_ep_sh.data,3);
    EOG.data.ext_blink_n_std(:,:,i_sub)     = std(temp.eog.EOG.EEG_ext_blink_ep_n.data,1,3);
    EOG.data.ext_blink_sh_std(:,:,i_sub)    = std(temp.eog.EOG.EEG_ext_blink_ep_sh.data,1,3);
    
    EOG.data.ext_resp_n_mean(:,:,i_sub)     = mean(temp.eog.EOG.EEG_ext_resp_ep_n.data,3);
    EOG.data.ext_resp_sh_mean(:,:,i_sub)    = mean(temp.eog.EOG.EEG_ext_resp_ep_sh.data,3);
    EOG.data.ext_resp_n_std(:,:,i_sub)      = std(temp.eog.EOG.EEG_ext_resp_ep_n.data,1,3);
    EOG.data.ext_resp_sh_std(:,:,i_sub)     = std(temp.eog.EOG.EEG_ext_resp_ep_sh.data,1,3);
    
    EOG.data.eeg_blink_n_mean(:,:,i_sub)    = mean(temp.eog.EOG.EEG_exp_blink_ep_n_ress.data,3);
    EOG.data.eeg_blink_sh_mean(:,:,i_sub)   = mean(temp.eog.EOG.EEG_exp_blink_ep_sh_ress.data,3);
    EOG.data.eeg_blink_n_std(:,:,i_sub)     = std(temp.eog.EOG.EEG_exp_blink_ep_n_ress.data,1,3);
    EOG.data.eeg_blink_sh_std(:,:,i_sub)    = std(temp.eog.EOG.EEG_exp_blink_ep_sh_ress.data,1,3);
    
    EOG.trialnum.ext_blink_n(i_sub)         = temp.eog.EOG.EEG_ext_blink_ep_n.trials;
    EOG.trialnum.ext_blink_sh(i_sub)        = temp.eog.EOG.EEG_ext_blink_ep_sh.trials;
    EOG.trialnum.ext_resp_n(i_sub)          = temp.eog.EOG.EEG_ext_resp_ep_n.trials;
    EOG.trialnum.ext_resp_sh(i_sub)         = temp.eog.EOG.EEG_ext_resp_ep_sh.trials;
    EOG.trialnum.eeg_blink_n(i_sub)         = temp.eog.EOG.EEG_exp_blink_ep_n_ress.trials;
    EOG.trialnum.eeg_blink_sh(i_sub)        = temp.eog.EOG.EEG_exp_blink_ep_sh_ress.trials;
    
    EOG.RESS_map(:,i_sub)                   = temp.eog.EOG.RESS.map(:,end);
    EOG.RESS_snr_ind(i_sub)                 = temp.eog.EOG.RESS.SNR_ind(end);
    EOG.RESS_snr_evo(i_sub)                 = temp.eog.EOG.RESS.SNR_evo(end);
    
    EOG.blinktimes_n{i_sub}                 = temp.eog.EOG.art.blink_trialtimes_n;
    EOG.blinktimes_sh{i_sub}                = temp.eog.EOG.art.blink_trialtimes_sh;
    
    EOG.blinkmat_n                          = [EOG.blinkmat_n temp.eog.EOG.art.blink_trialtimes_n];
    EOG.blinkmat_sh                         = [EOG.blinkmat_sh temp.eog.EOG.art.blink_trialtimes_sh];
    
    
    % TFA
    EOG.data.TFA_blink_noRESS_ind(:,:,:,i_sub) = temp.eog.EOG.TFA.data_induced;
    EOG.data.TFA_blink_noRESS_evo(:,:,:,i_sub) = temp.eog.EOG.TFA.data_evoked;
    EOG.data.TFA_blink_noRESS_ind_bc(:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.eog.EOG.TFA.data_induced, ...
        mean(temp.eog.EOG.TFA.data_induced(:,eeg_time2points(TFA.baseline(1),EOG.data.TFA_times):eeg_time2points(TFA.baseline(2),EOG.data.TFA_times),:,:),2)))-1);
    EOG.data.TFA_blink_noRESS_evo_bc(:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.eog.EOG.TFA.data_evoked, ...
        mean(temp.eog.EOG.TFA.data_evoked(:,eeg_time2points(TFA.baseline(1),EOG.data.TFA_times):eeg_time2points(TFA.baseline(2),EOG.data.TFA_times),:,:),2)))-1);
    EOG.data.TFA_blink_RESS_ind(1,:,:,i_sub) = temp.eog.EOG.TFA.data_RESS_induced;
    EOG.data.TFA_blink_RESS_evo(1,:,:,i_sub) = temp.eog.EOG.TFA.data_RESS_evoked;
    EOG.data.TFA_blink_RESS_ind_bc(1,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.eog.EOG.TFA.data_RESS_induced, ...
        mean(temp.eog.EOG.TFA.data_RESS_induced(:,eeg_time2points(TFA.baseline(1),EOG.data.TFA_times):eeg_time2points(TFA.baseline(2),EOG.data.TFA_times),:,:),2)))-1);
    EOG.data.TFA_blink_RESS_evo_bc(1,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.eog.EOG.TFA.data_RESS_evoked, ...
        mean(temp.eog.EOG.TFA.data_RESS_evoked(:,eeg_time2points(TFA.baseline(1),EOG.data.TFA_times):eeg_time2points(TFA.baseline(2),EOG.data.TFA_times),:,:),2)))-1);  
    clear temp
    
end

%% plotting blinktimes in trial
%figure; hist(EOG.blinkmat_n,100)
%figure; plot(sort(EOG.blinkmat_n))

pl.timewin = 100; % time window in ms
pl.data = nan(1,numel(EOG.time.ext_resp));
for i_t = 1:numel(EOG.time.ext_resp)
    tix = dsearchn(EOG.time.ext_resp',([-(pl.timewin/2) (pl.timewin/2)]+EOG.time.ext_resp(i_t))');
    pl.data(i_t)=sum((EOG.blinkmat_n >= EOG.time.ext_resp(tix(1))) & (EOG.blinkmat_n <= EOG.time.ext_resp(tix(2))));
end
figure; subplot(2,1,1)
plot(EOG.time.ext_resp,pl.data);
title(sprintf('counts of blinks | window = t +- %1.1f ms',pl.timewin/2));
xlabel('time in ms')
ylabel('count')

pl.data2 = abs(fft(detrend(pl.data),10000,2))*2/numel(pl.data);
pl.xdata= ((0:size(pl.data2,2)-1)/size(pl.data2,2)) * EOG.srate;

subplot(2,1,2)
plot(pl.xdata,pl.data2);
title(sprintf('frequency of counts of blinks | window = t +- %1.1f ms',pl.timewin/2));
xlabel('frequency in Hz')
ylabel('amplitude')
xlim([0 20])

%% do permutation test for blink distribution across experimental run
perm.timewin = [-3500 3500]; % time window for which blinks are anaysed
perm.bins = 60;
perm.n_reps = 10000; % number of repetitions
clear t

% first extract data across subjects
for i_sub = 1:numel(F.Subjects2Use)
    % only blinks in range are relevant
    t.data = EOG.blinktimes_n{i_sub}(EOG.blinktimes_n{i_sub}>=perm.timewin(1) & EOG.blinktimes_n{i_sub}<=perm.timewin(2));
    [t.hist.N,t.hist.edges,t.hist.bin] = histcounts(t.data,perm.bins,'BinLimits',perm.timewin);
    perm.bincenter = t.hist.edges(1:end-1)+(diff(t.hist.edges)./2);
    perm.o_data(i_sub,:)=t.hist.N;
end
% check plotting
% figure; bar(perm.bincenter,sum(perm.o_data))

% do permutation
perm.sh_data = nan(perm.n_reps,perm.bins);
for i_rep = 1:perm.n_reps
    % shuffle index for all subjects
    t.idx = cell2mat(cellfun(@(x) randperm(x), repmat(num2cell(perm.bins),1,numel(F.Subjects2Use)),'UniformOutput',false)');
    % create histogram according to this shuffling
    perm.sh_data(i_rep,:) = sum(perm.o_data(t.idx));
end

% find median and percentiles for each bin of empirical distribution
perm.sh_percentiles = nan(3,perm.bins);
for i_bin = 1:perm.bins
    perm.sh_percentiles(:,i_bin) = prctile(perm.sh_data(:,i_bin), [2.5 50 97.5]);
end

% plot results
pl.xlim = perm.timewin;
pl.li_col = [0.2 0.4 1];


clear figs
figs{1} =figure;
set(gcf,'Position',[100 100 600 200],'PaperPositionMode','auto')
pl.data = EOG.blinkmat_n; pl.data(pl.data<pl.xlim(1)|pl.data>pl.xlim(2)) = [];
histogram(pl.data,perm.bins,'BinLimits',pl.xlim, 'FaceColor', [0.2 0.4 1]);
%    title('distribution of events across experimental run')
hold on;
h.li = plot(perm.bincenter,perm.sh_percentiles);
set(h.li(1),'Color', pl.li_col, 'LineWidth',1.5, 'LineStyle',':')
set(h.li(2),'Color', pl.li_col, 'LineWidth',1.5, 'LineStyle','-')
set(h.li(3),'Color', pl.li_col, 'LineWidth',1.5, 'LineStyle',':')
xlim(pl.xlim)
ylim(get(gca,'YLim'))
xlabel('time in ms')
ylabel('count')
% title('Distirbution of blinks aligned to responses')

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'Blink_Disribution_by_response_withCI'};
for i_fig = 1:1
    print(figs{i_fig}, fullfile(sav.pathout,sav.filenames{i_fig}),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,sav.filenames{i_fig}),'fig')
    print(figs{i_fig},fullfile(sav.pathout,sav.filenames{i_fig}),'-depsc2', '-painters','-r300')
end

%% plotting blinktimes in trial
%figure; hist(EOG.blinkmat_n,100)
%figure; plot(sort(EOG.blinkmat_n))
pl.xlim = [-3500 3500];

clear figs
figs{1} =figure;
set(gcf,'Position',[100 100 600 200],'PaperPositionMode','auto')
pl.data = EOG.blinkmat_n; pl.data(pl.data<pl.xlim(1)|pl.data>pl.xlim(2)) = [];
histogram(pl.data,60,'BinLimits',pl.xlim, 'FaceColor', [0.2 0.4 1]);
%    title('distribution of events across experimental run')
xlim(pl.xlim)
ylim(get(gca,'YLim'))
xlabel('time in ms')
ylabel('count')
title('Distirbution of blinks aligned to responses')

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'Blink_Disribution_by_response'};
for i_fig = 1:1
    print(figs{i_fig}, fullfile(sav.pathout,sav.filenames{i_fig}),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,sav.filenames{i_fig}),'fig')
    print(figs{i_fig},fullfile(sav.pathout,sav.filenames{i_fig}),'-depsc2', '-painters','-r300')
end

%% plotting
pl.baseline = dsearchn(EOG.time.ext_blink',[-3500 -3000]');
pl.data = EOG.data.ext_blink_n_mean(1,:,:)-EOG.data.ext_blink_n_mean(2,:,:);
pl.data_bc = bsxfun(@minus, pl.data, mean(pl.data(:,pl.baseline(1):pl.baseline(2),:),2));
pl.data_ress = EOG.data.eeg_blink_sh_mean;

% filter
flt.fs = EOG.srate; %sampling rate
flt.tbw = 1; % transition bandwidth
flt.freqs =[-2*flt.tbw 2*flt.tbw]+(85/6);
flt.order = pop_firwsord('hamming', flt.fs, flt.tbw)*20;
[flt.b  flt.a] = firws(flt.order, flt.freqs / (flt.fs / 2), windows('hamming', flt.order + 1)); 
%freqz(flt.b,flt.a,size(pl.data,2),flt.fs);
% load surrogate data
EEG_s=EOG.surr_blink;
EEG_s=pop_select(EEG_s,'trial',1:numel(F.Subjects2Use));
EEG_s.data = pl.data_ress;
% EEG_s_f = firfilt(EEG_s, flt.b); % firfilt
EEG_s_f2 = pop_eegfiltnew(EEG_s, flt.freqs(1),flt.freqs(2), flt.order, 0, [], 0); % filter nr 2
% pop_eegplot(EEG_s_f2,1,1,1)
% pop_eegplot(EEG_s,1,1,1)

EEG_s_gab = eegF_Gabor(EEG_s,85/6,diff([-2*flt.tbw 2*flt.tbw]));

pl.data_ress_filt = EEG_s_f2.data;
pl.data_ress_filt_hilb = pl.data_ress_filt;
for i_sub = 1:numel(F.Subjects2Use)
    pl.data_ress_filt_hilb(1,:,i_sub) = abs(hilbert(pl.data_ress_filt(1,:,i_sub)));
end
%figure; plot(EEG_s.times,EEG_s.data(1,:,1));hold on;  plot(EEG_s.times,pl.data_ress_filt(1,:,1));plot(EEG_s.times,pl.data_ress_filt_hilb(1,:,1))
pl.baseline = dsearchn(EOG.time.ext_blink',[-3000 -2500]');
pl.data_ress_filt_hilb_bc = bsxfun(@minus, pl.data_ress_filt_hilb, mean(pl.data_ress_filt_hilb(:,pl.baseline(1):pl.baseline(2),:),2));
pl.data_ress_gab = EEG_s_gab.data;
pl.data_ress_gab_bc = bsxfun(@minus, pl.data_ress_gab, mean(pl.data_ress_gab(:,pl.baseline(1):pl.baseline(2),:),2));
pl.data_ress_bc = bsxfun(@minus, pl.data_ress, mean(pl.data_ress(:,pl.baseline(1):pl.baseline(2),:),2));
pl.data_ress_filt_bc = bsxfun(@minus, pl.data_ress_filt, mean(pl.data_ress_filt(:,pl.baseline(1):pl.baseline(2),:),2));

% plotting
clear figs
figs{1} = figure;
set(gcf,'Position',[100 100 600 200],'PaperPositionMode','auto')
pl.ydata1 = mean(pl.data,3);
% pl.ydata2 = [mean(pl.data_ress_filt_hilb_bc,3); mean(pl.data_ress_gab_bc,3); mean(pl.data_ress_bc,3); mean(pl.data_ress_filt_bc,3)];
pl.ydata2 = [mean(pl.data_ress_filt_hilb,3); mean(pl.data_ress_gab,3); mean(pl.data_ress,3); mean(pl.data_ress_filt,3)];
pl.ydata2 = [mean(pl.data_ress_gab,3); mean(pl.data_ress,3)];

% pl.ydata1 = mean(pl.data(:,:,1),3);
% pl.ydata2 = [mean(pl.data_ress_filt_hilb(:,:,1),3); mean(pl.data_ress_gab(:,:,1),3); mean(pl.data_ress(:,:,1),3); mean(pl.data_ress_filt(:,:,1),3)];

% [h.AX,h.line1,h.line2] = plotyy(EEG_s.times,pl.ydata1,repmat(EEG_s.times,size(pl.ydata2,1),1)',pl.ydata2','LineWidth',2);
colororder([ 0 0.5 0.8; 0.8 0 0 ])
yyaxis left
h.line1 = plot(EEG_s.times,pl.ydata1,'LineWidth',2);
grid on
set(gca,'ylim',[-1 1]*max(abs(get(gca,'ylim'))),'YTickMode','auto','xlim',[-3500 3500])
ylabel('amplitude of blinks in \muV');
set(h.line1(1),'Color',[0 0.5 0.8])

yyaxis right
h.line2 = plot(EEG_s.times,pl.ydata2);
t.max = max(max(max(abs([pl.ydata2]))));
%t.sc = ceil(t.max*(10^(fix(1/t.max)-1)))/10^(fix(1/t.max)-1);
if t.max<1
    t.sc = ceil(t.max*10^numel(num2str(fix(1/t.max))))/(10^numel(num2str(fix(1/t.max))));
else
    t.sc = ceil(t.max*1)/1;
end
set(gca,'ylim',[-1 1]*t.sc,'YTickMode','auto')
ylabel('amplitude of SSVEP in \muV');
set(h.line2(1),'Color',[0.8 0 0])
set(h.line2(2),'Color',[1 .75 .75],'LineStyle','-')

xlabel('time in ms')
% axes(h.AX(1)); ylabel('blinks'' amplitude in \muV');
% axes(h.AX(2)); ylabel('SSVEP amplitude modulation in \muV');
% legend([h.line1; h.line2],{'blink';'hilbert';'gabor';'RESS signal';'RESS signal filt'})
% axes(h.AX(2)); ylabel('SSVEP amplitude in \muV');
legend([h.line1; h.line2],{'blink';'SSVEP amplitude';'raw RESS signal'},'Location','SouthEast')


sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'Blink_SSVEP_Amp_TimeCourse'};
for i_fig = 1:1
    print(figs{i_fig}, fullfile(sav.pathout,sav.filenames{i_fig}),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,sav.filenames{i_fig}),'fig')
    print(figs{i_fig},fullfile(sav.pathout,sav.filenames{i_fig}),'-depsc2', '-painters','-r300')
end

%% actual plotting no RESS data blink
% figure; topoplot([],TFA.electrodes, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);


% plotting parameters
% pl.elec2plot = {'Oz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; % vis alpha
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
% pl.elec2plot = {'POz'}; % vis alpha II
% pl.elec2plot = {'C3';'CP3'}; % motor alpha/beta
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; sav.chan_add = 'MotorLarge';% motor alpha/beta II
% pl.elec2plot = {'F1';'F2';'Fz'}; sav.chan_add = 'Frontal'; % frontal
% pl.elec2plot = {'FP1';'FP2'}; % frontal
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; sav.chan_add = 'Parietal'; % parietal
% pl.elec2plot = {'F5';'F3';'F1';'FC3'}; sav.chan_add = 'PreMotor'; % pre-motor
% pl.elec2plot = {'FC3';'FC3';'FC1';'FCz';'Fz';'F5';'F3';'F1';'AF7';'AF3';'AFz';'FP1';'FPz'}; sav.chan_add = 'CentroFrontalLeft'; % pre-motor
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({EOG.data.TFA_electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% pl.xlims=[-2500 2500]; % index time 2 plot
pl.xlims=[-3500 3500]; % index time 2 plot
[t.t t.ind1]=min(abs(EOG.data.TFA_times-pl.xlims(1)));
[t.t t.ind2]=min(abs(EOG.data.TFA_times-pl.xlims(2)));

% pl.flims = EOG.data.TFA_freqs([1 end]); % index frequency 2 plot
% pl.flims = [4 35];
% pl.flims = [8 44];
pl.flims = [8 EOG.data.TFA_freqs([end])]; % index frequency 2 plot
[t.t t.find1]=min(abs(EOG.data.TFA_freqs-pl.flims(1)));
[t.t t.find2]=min(abs(EOG.data.TFA_freqs-pl.flims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
clear figs fig


%%%%%%%%%%%%%
% raw induced no RESS
pl.data_ind = squeeze(nanmean(nanmean(EOG.data.TFA_blink_noRESS_ind(:,:,pl.elec2plot_i,pl.subs2use),3),4));
pl.data_evo = squeeze(nanmean(nanmean(EOG.data.TFA_blink_noRESS_evo(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));
figs{1}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
subplot(2,1,1)
imagesc(EOG.data.TFA_times,EOG.data.TFA_freqs,pl.data_ind,t.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for induced activity\n for channel [%s]', vararg2str(pl.elec2plot)), 'FontSize',8)
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
set(gca,'FontSize',8)

subplot(2,1,2)
imagesc(EOG.data.TFA_times,EOG.data.TFA_freqs,pl.data_evo,t.clims2)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for evoked activity\n for channel [%s]', vararg2str(pl.elec2plot)), 'FontSize',8)
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
set(gca,'FontSize',8)
% colormap('jet')

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),EOG.data.TFA_electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});



% normalized induced no RESS
pl.data_ind = squeeze(nanmean(nanmean(EOG.data.TFA_blink_noRESS_ind_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
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
imagesc(EOG.data.TFA_times,EOG.data.TFA_freqs,pl.data_ind,t.clims1)
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

h.sp2=subplot(2,1,2);
for i_freq = 1:numel(EOG.data.TFA_freqs)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(nanmean(EOG.data.TFA_blink_noRESS_ind_bc(i_freq,:,pl.elec2plot_i,pl.subs2use),3))');
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
imagesc(EOG.data.TFA_times,EOG.data.TFA_freqs,abs(log10(TFA.pvals_induced_bc)),t.clim)
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
topoplot(find(pl.elec2plot_i),EOG.data.TFA_electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});


% normalized evoked noRESS
pl.data_evo = squeeze(nanmean(nanmean(EOG.data.TFA_blink_noRESS_evo_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims2=[-1 1]*max(max(abs(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))));
figs{3}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
h.sp1=subplot(2,1,1);
imagesc(EOG.data.TFA_times,EOG.data.TFA_freqs,pl.data_evo,t.clims2)
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
for i_freq = 1:numel(EOG.data.TFA_freqs)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(nanmean(EOG.data.TFA_blink_noRESS_evo_bc(i_freq,:,pl.elec2plot_i,pl.subs2use),3))');
    TFA.pvals_evoked_bc(i_freq,:)=tt.p;
end
t.data = abs(log10(TFA.pvals_evoked_bc(t.find1:t.find2,:)));
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
imagesc(EOG.data.TFA_times,EOG.data.TFA_freqs,abs(log10(TFA.pvals_evoked_bc)),t.clim)
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
topoplot(find(pl.elec2plot_i),EOG.data.TFA_electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'Blink_TFA_Amp_RAW_EvoIndu';'Blink_TFA_Amp_BC_Indu';'Blink_TFA_Amp_BC_Evo'};
for i_fig = 1:3
    print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'fig')
    print(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-depsc2', '-painters','-r300')
end


%% plot all lines of interest into one graphics (bounded line)
clear figs h
pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'RESS_evo_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
    [14.16667 14.16667] {'RESS'} 'RESS_ind_bc' [0.5 0.5 0.5] 2 sprintf('vis ind\nSSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'noRESS_ind_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [17 20] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'noRESS_ind_bc' [0 0 0.8] 2 'vis beta' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'noRESS_ind_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    };

pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'RESS_evo_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
    [14.16667 14.16667] {'RESS'} 'RESS_ind_bc' [0.5 0.5 0.5] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'noRESS_ind_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [17 20] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'noRESS_ind_bc' [0 0 0.8] 2 'vis beta' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'noRESS_ind_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    };



pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
% 
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % subjects discarded due to low trial number version 2


% pl.xlim = [-3250 3250];
pl.xlim = [-4000 4000];
pl.xlim = [-3500 3500];
% pl.xlim = [-3000 3000];
pl.permut_n = 10000;

pl.toponum = [ceil(sqrt(size(pl.parameters,1))) round(sqrt(size(pl.parameters,1)))];
pl.plpos = repmat([1:pl.toponum(2)*3],pl.toponum(1),1)+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2)*3);
pl.topopos = repmat([pl.toponum(2)*3+1:pl.toponum(2)*4],pl.toponum(1),1)'+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2))';

h.pl = [];
figs{1}=figure;
set(gcf,'Position',[100 100 800 350],'PaperPositionMode','auto')
h.sp(1)=subplot(pl.toponum(1),pl.toponum(2)*4,pl.plpos(:));
for i_pl = 1:size(pl.parameters,1)
    % index frequencies
    [t.t t.ind1]=min(abs(EOG.data.TFA_freqs-pl.parameters{i_pl,1}(1)));
    [t.t t.ind2]=min(abs(EOG.data.TFA_freqs-pl.parameters{i_pl,1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(EOG.data.TFA_times-pl.xlim(1)));
    [t.t t.ind4]=min(abs(EOG.data.TFA_times-pl.xlim(2)));
    
    % index electrodes
    switch  pl.parameters{i_pl,7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({EOG.data.TFA_electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data=squeeze(nanmean(nanmean(EOG.data.TFA_blink_%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.elec2plot_i,pl.subs2use),1),3));',pl.parameters{i_pl,3});
        case 'RESS'
            com=sprintf('pl.data=squeeze(nanmean(EOG.data.TFA_blink_%s(1,t.ind1:t.ind2,t.ind3:t.ind4,pl.subs2use),2));',pl.parameters{i_pl,3});
        
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
    [h.l(i_pl), h.p(i_pl)] = boundedline(EOG.data.TFA_times(t.ind3:t.ind4),pl.mdata,pl.ddata(1,:));
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
set(gca,'ylim',[-(1+0.15+0.05*size(pl.parameters,1)) 1.1].*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
pl.maxdata = max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false)));
for i_pl = 1:size(pl.parameters,1)
    pl.y = -1*pl.maxdata*(1+0.05+0.05*i_pl);
    
%     plot(TFA.time(t.ind3:t.ind4),repmat(pl.y,1,numel(TFA.time(t.ind3:t.ind4))),'Color','k','LineWidth',0.1)
    plot(EOG.data.TFA_times(t.ind3:t.ind4),pl.tdata2(i_pl,:).*pl.y,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',2)
    
end



legend(h.l, pl.parameters(:,6),'location','NorthWest','FontSize', 8)
grid on
box on
hline(0,'k')
xlabel('time in ms')
ylabel('amplitude modulation in %')

for i_pl = 1:size(pl.parameters,1)
    subplot(pl.toponum(1),pl.toponum(2)*4,pl.topopos(i_pl))
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({EOG.data.TFA_electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
%     topoplot([],TFA.electrodes(1:64),'style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',2});
    topoplot(find(pl.elec2plot_i),EOG.data.TFA_electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
        'emarker2',{find(pl.elec2plot_i),'o','r',3,1});
    title(sprintf('%s\n[%1.0f %1.0f]Hz',pl.parameters{i_pl,6},pl.parameters{i_pl,1}),...
        'Color',pl.parameters{i_pl,4})
end

% sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
% sav.filenames = {'Blink_AllSignals_Amp_Timecourse'};
% for i_fig = 1:1
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s,',sav.filenames{i_fig})),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s,',sav.filenames{i_fig})),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s,',sav.filenames{i_fig})),'-depsc2', '-painters','-r300')
% end
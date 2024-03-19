clear all;
%% script in order to determine ERD, see ERD pattern
% cuts data into pre and post stimulus epochs
% calculates pre, post and difference fft spectra
% additionally calculates time frequency analysis and wavelet sepctrum
% implemented for different references and elecrodes to plot


%% %%%%%%%%%%%%%%%definitions%%%%%%%%%%% %%
F.Subjects =        ['10']; % VP 05 % VP 07 (nice EMG) % VP 08 (nice Beta topo) % VP 09 (nice alpha) % VP 10 (best so far?)
F.EpochRangeEMG =   [-3 3];     % Timerange of epoch in s
F.EpochRangeTFA =   [-3 3];
F.Trigger =         {'31'};
F.Chans2Plot =      {'C3';'POz'};   % chans to plot
F.Chans2Plot =      {{'C3';'CP3';'C1';'CP1'};{'POz'}};   % chans to plot  
F.FiltFreq =        [];
F.ScalePlot =       [3 25 -1 7];
F.RefChan =         {'avg'}; % 'noreref' = no rereference
%F.RefChan =         {'C3';'C4'}; % 'noreref' = no rereference
F.TFAFreqs =        [2:0.1:30];
F.TFABaseline =     [-2 -1.5]; % [-2.25 -2];
F.MuRange =         [8 14];
F.EEGChans =        [1:64];
F.EMGChans =        71:72;
F.EMGBaseline =     [-1.75 -1.5]; % in s [-3 -2.5]


F.RAWPath =         ['F:\work\data\SSVEP_volmov\EEG\RAW\'];
F.ChanlocFile =     'C:\Users\HP-User\matlab\Skripte\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.FileName =        ['_motor'];
F.Prefix =          'VP';

Disp.xlim =         [-2000 2000];

%% read-in data
fprintf(1,'reading in data %s%s%s%s.bdf\n', F.RAWPath,F.Prefix,F.Subjects,F.FileName)
EEG=pop_biosig(sprintf('%s%s%s%s.bdf', F.RAWPath,F.Prefix,F.Subjects,F.FileName));
% pop_eegplot(EEG,1,1,1);
% chanlocs
EEG.setname=[F.Prefix '_' F.Subjects 'alpha'];
% pop_eegplot(EEG,1,1,1)
% change locations
EEG.chanlocs=pop_chanedit(EEG.chanlocs,'load',{F.ChanlocFile,'filetype','besa (elp)'}); % Load Channel Locations
% EEG=pop_chanedit(EEG, 'lookup','C:\Users\HP-User\matlab\software\toolboxes\eeglab13_3_2b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp');


EMG = pop_select( EEG,'channel',F.EMGChans);
EEG = pop_select( EEG,'channel',F.EEGChans);
% pop_eegplot(EEG,1,1,1)
%figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);

%% filter
if ~isempty(F.FiltFreq)
    Wn2=F.FiltFreq/(EEG.srate/2);
    [b2,a2]=butter(3,Wn2,'stop');
    for i_el = 1:size(EEG.data,1)
        EEG.data(i_el,:) = filtfilt(b2,a2,double(EEG.data(i_el,:)));
    end
end


%% reref
if ~any(strcmpi(F.RefChan{1},'noreref'))
    % determine reference channels
    if ~any(strcmpi(F.RefChan{1},'avg'))
        Temp.RefChan = [];
        for i_el = 1:length(F.RefChan{1})
            Temp.RefChan(i_el) = find(strcmp({EEG.chanlocs(1,:).labels},F.RefChan{i_el}));
        end
        EEG_t1 = pop_reref( EEG, Temp.RefChan);
    else
        EEG_t1 = pop_reref( EEG, []); % average
    end
else
    EEG_t1 = EEG;
    % pop_eegplot(EEG_t1,1,1,1)
end

%% epoch
EEG_tfa = pop_epoch( EEG_t1, F.Trigger, F.EpochRangeTFA, 'epochinfo', 'yes');
EEG_tfa = pop_rmbase( EEG_tfa, []);
EEG_tfa = eeg_checkset( EEG_tfa );

EMG_ep = pop_epoch( EMG, F.Trigger, F.EpochRangeEMG, 'epochinfo', 'yes');
EMG_ep = pop_rmbase(EMG_ep, F.EMGBaseline);
EMG_ep = eeg_checkset( EMG_ep );

% filter
EEG_behav_f = pop_eegfiltnew(EEG_t1, 0.25, 8, EEG_t1.srate*8, 0, [], 0);
EEG_ep = pop_epoch( EEG_behav_f, F.Trigger, F.EpochRangeEMG, 'epochinfo', 'yes');
EEG_ep = pop_rmbase(EEG_ep, F.EMGBaseline*1000);
EEG_ep = eeg_checkset( EEG_ep );
% pop_eegplot(EEG_ep,1,1,1)

%% preallocate memory

% TFAData=nan(numel(F.Chans2Plot),length(F.TFAFreqs),size(EEG_tfa.data,2));
% % TFAData =     1. Dim: channels
% %               2. Dim: frequencies
% %               3. Dim: datapoints

TFA.Data=nan(numel(F.EEGChans),length(F.TFAFreqs),size(EEG_tfa.data,2));
% TFA.Data =     1. Dim: channels
%               2. Dim: frequencies
%               3. Dim: datapoints



%% calculate TFA
%fprintf(1,'calculating TFA at %1.0f channels',numel(F.Chans2Plot))
fprintf(1,'calculating TFA at %1.0f channels ',numel(F.EEGChans))
for i_el=1:numel(F.EEGChans) %loop electrodes
    %Temp.Chan = find(strcmp({EEG_tfa.chanlocs(:).labels},F.Chans2Plot{i_el}));
    %fprintf(1,' %s',F.Chans2Plot{i_el})
    %[TFA.Data(i_el,:,:),TFA.t,TFA.f]=traces2TFA_4(squeeze(EEG_tfa.data(Temp.Chan,:,:)),F.TFAFreqs,EEG_tfa.srate,7);
    %fprintf(1,' %s',EEG_tfa.chanlocs(i_el).labels)
    fprintf(1,'%1.0f',mod(i_el,10))
    [TFA.Data(i_el,:,:),TFA.t,TFA.f]=traces2TFA_4(squeeze(EEG_tfa.data(i_el,:,:)),F.TFAFreqs,EEG_tfa.srate,7);
end
fprintf(1,'...done\n')
        
    

%% display TFA-results
% loop for electrodes 2 plot
pl.el2plot={};
for i_elc = 1:numel(F.Chans2Plot)
    pl.el2plot{i_elc}=cellfun(@(x) find(strcmpi({EEG.chanlocs.labels},x)), F.Chans2Plot{i_elc})';
end

TFA.Data_c=TFA.Data-repmat(...
    mean(TFA.Data(:,:,eeg_time2points(F.TFABaseline(1)*1000,EEG_tfa.times):eeg_time2points(F.TFABaseline(2)*1000,EEG_tfa.times)),3),...
    [1,1,size(TFA.Data,3)]);
for i_chans2plot = 1:length(F.Chans2Plot)
    fig1 = figure;
    set(fig1,'WindowStyle','docked')
    [~, TEMP.xlim1] = min(abs(EEG_tfa.times-Disp.xlim(1)));
    [~, TEMP.xlim2] = min(abs(EEG_tfa.times-Disp.xlim(2)));
    % loop for references
    
    
    % plot tfa
    figure(fig1)
    subplot(2,1,1)
    % clims = [-1 1]*max(max(abs(squeeze(TFA.Data(i_chans2plot,:,:,i_ref2plot)))));
    clims = [min(min(squeeze(mean(TFA.Data(pl.el2plot{i_chans2plot},:,TEMP.xlim1:TEMP.xlim2),1)))) ...
        max(max(squeeze(mean(TFA.Data(pl.el2plot{i_chans2plot},:,TEMP.xlim1:TEMP.xlim2),1))))];
    imagesc(EEG_tfa.times,TFA.f,squeeze(mean(TFA.Data(pl.el2plot{i_chans2plot},:,:),1)),clims)
    set(gca,'YDir','normal')
    title(sprintf('tfa, for channel %s referenced to %s',...
        vararg2str(F.Chans2Plot{i_chans2plot}),vararg2str(F.RefChan{:})))
    colorbar
    xlabel('time in ms')
    ylabel('frequency in Hz')
    xlim(Disp.xlim)
    
    % plot tfa baseline corrected
    subplot(2,1,2)
    clims = [-1 1]*max(max(abs(mean(TFA.Data_c(pl.el2plot{i_chans2plot},:,TEMP.xlim1:TEMP.xlim2),1))));
    % clims = [-1 1]*abs(min(min(TFA2plotB(:,TEMP.xlim1:TEMP.xlim2))));
    % clims = [min(min(TFA2plotB)) max(max(TFA2plotB))];
    imagesc(EEG_tfa.times,TFA.f,squeeze(mean(TFA.Data_c(pl.el2plot{i_chans2plot},:,:),1)),clims)
    set(gca,'YDir','normal')
    title(sprintf('tfa, baseline corrected - for channel %s referenced to %s; baseline [%1.0f %1.0f]ms',...
        vararg2str(F.Chans2Plot{i_chans2plot}),vararg2str(F.RefChan{:}),F.TFABaseline*1000 ))
    colorbar
    xlabel('time in ms')
    ylabel('frequency in Hz')
    xlim(Disp.xlim)
    
end

%% display TFA topo
pl.timerange = [-500 1000];
pl.freqrange = [8 10];
% pl.freqrange = [15 23];
[t.t pl.tind1]=min(abs(EEG_tfa.times-pl.timerange(1)));
[t.t pl.tind2]=min(abs(EEG_tfa.times-pl.timerange(2)));
[t.t pl.find1]=min(abs(TFA.f-pl.freqrange(1)));
[t.t pl.find2]=min(abs(TFA.f-pl.freqrange(2)));

pl.data1 = TFA.Data;
pl.data2 = TFA.Data_c;

pl.topodata1=mean(mean(pl.data1(:,pl.find1:pl.find2,pl.tind1:pl.tind2),2),3);
pl.topodata2=mean(mean(pl.data2(:,pl.find1:pl.find2,pl.tind1:pl.tind2),2),3);

topofig = figure;
set(topofig,'WindowStyle','docked')

subplot(1,2,1)
topoplot( pl.topodata1, EEG_tfa.chanlocs(F.EEGChans), 'chaninfo', EEG_tfa.chaninfo, ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',[0 1]*max(abs(pl.topodata1)));
title(sprintf('raw amplitde t[%1.0f %1.0f]ms;\nf[%1.1f %1.1f]Hz',pl.timerange,pl.freqrange))
colorbar

subplot(1,2,2)
topoplot( pl.topodata2, EEG_tfa.chanlocs(F.EEGChans), 'chaninfo', EEG_tfa.chaninfo, ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',[-1 1]*max(abs(pl.topodata2)));
title(sprintf('baseline corrected amplitde t[%1.0f %1.0f]ms;\nf[%1.1f %1.1f]Hz',pl.timerange,pl.freqrange))
colorbar

%% display EMG
pl.data = mean(EMG_ep.data(2,:,:)-EMG_ep.data(1,:,:),3);

figure; plot(EMG_ep.times,pl.data);
title('EMG response')
xlabel('time in ms')
ylabel('amplitude in \muV')
grid on

%% display EEG
figure; 
pop_plottopo(EEG_ep, 1:64 , 'BDF file epochs', 0, 'ydir',-1);

figure; 
EEG_ept = pop_select( EEG_ep,'time',[-0.2 0.5] );
pop_plottopo(EEG_ept, 1:64 , 'BDF file epochs', 0, 'ydir',-1);


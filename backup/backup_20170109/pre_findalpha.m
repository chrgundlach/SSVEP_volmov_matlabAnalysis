clear all;
%% script in order to determine ERD
% cuts data into pre and post stimulus epochs
% calculates pre, post and difference fft spectra
% additionally calculates time frequency analysis and wavelet sepctrum
% implemented for different references and elecrodes to plot

%% %%%%%%%%%%%%%%%definitions%%%%%%%%%%% %%
F.Subjects =        ['01'];
F.EpochRangeFFT =   [0 2];     % Timerange of epoch in s
F.Trigger =         {'26';'25'};
F.Trigger_con =     {'eyes closed';'eyes open'};
F.Blocks =          3;
F.TrigPerBlock =    [30];
F.Chans2Plot =      {'C3';'Oz';'O2'};   % chans to plot            
F.FFTRes =          25000;
F.FiltFreq =        [];
F.ScalePlot =       [3 25 -1 7];
F.RefChan =         {'avg'}; % 'noreref' = no rereference
%F.RefChan =         {'C3';'C4'}; % 'noreref' = no rereference
F.TFAFreqs =        [5:0.1:30];
F.TFASpectrum =     [0 2];
F.AlphaRange =      [8 14];
F.EEGChans =        [1:64];

F.RAWPath =         ['F:\work\data\SSVEP_volmov\EEG\RAW\'];
F.ChanlocFile =     'C:\Users\HP-User\matlab\Skripte\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.FileName =        ['_alpha'];
F.Prefix =          'VP';

Disp.xlim =         [0 2000];

%% read-in data
fprintf(1,'reading in data %s%s%s%s.bdf\n', F.RAWPath,F.Prefix,F.Subjects,F.FileName)
EEG=pop_biosig(sprintf('%s%s%s%s.bdf', F.RAWPath,F.Prefix,F.Subjects,F.FileName));
% pop_eegplot(EEG,1,1,1);
% chanlocs
EEG.setname=[F.Prefix '_' F.Subjects 'alpha'];
% pop_eegplot(EEG,1,1,1)
% change locations
EEG.chanlocs=pop_chanedit(EEG.chanlocs,'load',{F.ChanlocFile,'filetype','besa (elp)'}); % Load Channel Locations
EEG = pop_select( EEG,'channel',F.EEGChans);

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
        %EEG_t1 = pop_reref( EEG, [],'exclude',[F.EEGChans+1:EEG.nbchan] ); % average
        EEG_t1 = pop_reref( EEG, []); % average
    end
else
    EEG_t1 = EEG;
    % pop_eegplot(EEG_t1,1,1,1)
end

%% epoch for two triggers
for i_ep = 1:numel(F.Trigger)
    EEG_fft{i_ep} = pop_epoch( EEG_t1, F.Trigger(i_ep), F.EpochRangeFFT, 'epochinfo', 'yes');
    EEG_fft{i_ep} = pop_select( EEG_fft{i_ep},'trial',setdiff(1:EEG_fft{i_ep}.trials ,1:F.TrigPerBlock:F.Blocks*F.TrigPerBlock) );
    EEG_fft{i_ep} = pop_rmbase( EEG_fft{i_ep}, []);
    EEG_fft{i_ep} = eeg_checkset( EEG_fft{i_ep} );
end

%% preallocate memory
FFTData = nan(2,size(EEG_fft{1}.data,3),numel(F.EEGChans), F.FFTRes);
% FFTData =     1. Dim: eyes open /eyes closed
%               2. Dim: epoch
%               3. Dim: Channels
%               4. Dim: Datapoints
%               5. Dim: references


%% calculate FFT
%determine channels to plot

% detrend first
EEG_fft{1}= eegF_Detrend(EEG_fft{1});
EEG_fft{2}= eegF_Detrend(EEG_fft{2});
% loop for each epoch
for i_epoch = 1:size(EEG_fft{1}.data,3)
    FFTData(1,i_epoch,:,:) = abs(fft(EEG_fft{1}.data(F.EEGChans,:,i_epoch),F.FFTRes,2))*2/size(EEG_fft{1}.data,2);
end
for i_epoch = 1:size(EEG_fft{2}.data,3)
    FFTData(2,i_epoch,:,:) = abs(fft(EEG_fft{2}.data(F.EEGChans,:,i_epoch),F.FFTRes,2))*2/size(EEG_fft{2}.data,2);
end
% pop_eegplot(EEG_fft{1},1,1,1)
% pop_eegplot(EEG_fft{2},1,1,1)

%% create subplot indices
Subp.VerNum = 3;
Subp.HorNum = length(F.RefChan);

%% calculate of FFT means
% mean for epochs
FFTData = nanmean(FFTData,2);

% calculate difference of pre and post
FFTData(3,:,:,:,:) = ([FFTData(1,:,:,:,:) - FFTData(2,:,:,:,:)]);

%% find maxima
% calculate scale
xScale = ((0:size(FFTData,4)-1)/size(FFTData,4)) * EEG.srate;
[~, SPRange1]=min(abs(xScale-F.AlphaRange(1)));
[~, SPRange2]=min(abs(xScale-F.AlphaRange(2)));

SPRange = (SPRange1 : SPRange2);
pl.el2plot=cellfun(@(x) find(strcmpi({EEG.chanlocs.labels},x)), F.Chans2Plot)';

FFTAlpha = zeros(length(F.Chans2Plot),1);
for i_el = 1:numel(pl.el2plot)
    
    [I Freq] = max(squeeze(FFTData(3,1,pl.el2plot(i_el),SPRange)));
    FFTAlpha(i_el) = xScale(SPRange(Freq));
    
end


%% plot fft data
pl.el2plot=cellfun(@(x) find(strcmpi({EEG.chanlocs.labels},x)), F.Chans2Plot)';
for i_el = 1:length(F.Chans2Plot)
    fftfig = figure;
    set(fftfig,'WindowStyle','docked')
    %  red = pre, blue = post
    plot(xScale,zeros(1,size(FFTData,4)), '--k');
    hold on;
    hpl = plot(xScale,squeeze(FFTData(1,1,pl.el2plot(i_el),:))','r', ...
        xScale,squeeze(FFTData(2,1,pl.el2plot(i_el),:))', 'b', ...
        xScale,squeeze(FFTData(3,1,pl.el2plot(i_el),:))', 'g', 'LineWidth',2);
    % axis(F.ScalePlot);
    xlim([F.ScalePlot(1) F.ScalePlot(2)])
    set(gca, 'Ylim',[F.ScalePlot(3) max(get(gca, 'Ylim'))])
    title({'average FFT-spectra for pre- and post stimulus time window';['one participant; electrode ' F.Chans2Plot{i_el} ' referenced to ' vararg2str(F.RefChan{:})]});
    
    xlabel('frequency (Hz)');
    ylabel('amplitude in \muV');
    vline(FFTAlpha(i_el),'k', ['individual alpha frequency '  num2str(FFTAlpha(i_el)) ' Hz'])
    legend(hpl,F.Trigger_con{1},F.Trigger_con{2}, 'diff');
    
end

%% plot topo
% index frequencies
pl.freqrange = [-0.5 0.5];
pl.freq2plot=3;
[t.t t.ind1] = min(abs(xScale-(FFTAlpha(pl.freq2plot)+pl.freqrange(1))));
[t.t t.ind2] = min(abs(xScale-(FFTAlpha(pl.freq2plot)+pl.freqrange(2))));

pl.data=mean(FFTData(:,:,:,t.ind1:t.ind2),4);
pl.con=[F.Trigger_con; {'diff'}];

topofig = figure;
set(topofig,'WindowStyle','docked')
t.clim=[0 max(max(pl.data(1:2,:,:))); 0 max(max(pl.data(1:2,:,:))); [-1 1]*max(abs(pl.data(3,:,:)))];
for i_pl=1:3
    subplot(1,3,i_pl)
    
    topoplot( squeeze(pl.data(i_pl,:,:)), EEG.chanlocs(F.EEGChans), 'chaninfo', EEG.chaninfo, ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',t.clim(i_pl,:));
    title(sprintf('%s; freq %1.2fHz + [%1.2f %1.2f]Hz',pl.con{i_pl},FFTAlpha(pl.freq2plot),pl.freqrange))
    colorbar
end


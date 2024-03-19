%% first analysis steps

clear all
%% load file

%% parameters
F.WorkingPath           = 'F:\work\data\SSVEP_volmov\EEG\';
F.File                  = 'VP02_exp';
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.ChanlocFile           = 'C:\Users\HP-User\matlab\Skripte\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.RespTrigger           = 30;
F.RespEpoch             = [-2.5 2.5]; % in s
F.SSVEPTrigger          = {'15'};
F.SSVEPEpoch            = [0 6]; % in s

F.ExpEpochTrigger       = {'11','12'};

%% preprocessing
%load file
EEG=pop_biosig(sprintf('%sraw\\%s.bdf', F.WorkingPath,F.File));
% pop_eegplot(EEG,1,1,1);

% chanlocs
EEG.chanlocs=pop_chanedit(EEG.chanlocs,'load',{F.ChanlocFile,'filetype','besa (elp)'}); % Load Channel Locations
%figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);

% rereference
EEG = pop_reref( EEG, [],'exclude',[F.EEGChans(end)+1:EEG.nbchan] ); % average
% EEG = pop_reref( EEG, 69); % earlobe

% pop_eegplot(EEG,1,1,1,1)

%% behavior
EEG_behav = pop_epoch( EEG, {F.RespTrigger}, F.RespEpoch,'epochinfo', 'yes');
EEG_behav = pop_rmbase( EEG_behav, [-2250 -1750]);

% filter
EEG_behav_f = pop_eegfiltnew(EEG_behav, 0.25, 8, EEG.srate*diff(F.RespEpoch)*2, 0, [], 0);    

% plot EEG
pl.elec2plot={'C3';'CP3'};
pl.elec2plot={'Cz'};
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({EEG_behav_f.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
figure;
plot(EEG_behav_f.times,squeeze(nanmean(nanmean(EEG_behav_f.data(pl.elec2plot_i,:,:),1),3)))
title(sprintf('average ERP - electrodes %s',vararg2str(pl.elec2plot)))
xlabel('time in ms')
ylabel('amplitude in \muV')
set(gca,'ydir','reverse')

% plot EEG in scalp
figure; pop_plottopo(EEG_behav_f, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);

% plot EMG
figure;
plot(EEG_behav.times,squeeze(nanmean(EEG_behav.data(F.EMGChans(1),:,:)-EEG_behav.data(F.EMGChans(2),:,:),3)));
title(sprintf('average EMG - electrodes %s'))
xlabel('time in ms')
ylabel('amplitude in \muV')
set(gca,'ydir','reverse')


%% SSVEP
EEG_SSVEP = pop_epoch( EEG, F.SSVEPTrigger, F.SSVEPEpoch,'epochinfo', 'yes');
EEG_SSVEP = pop_rmbase( EEG_SSVEP, []);
EEG_SSVEP = eegF_Detrend(EEG_SSVEP);

% discard erroneous trials
t.tr2discard = unique(cell2mat({EEG_SSVEP.event(cell2mat({EEG_SSVEP.event.type})== 11 | cell2mat({EEG_SSVEP.event.type})== 12).epoch}));
EEG_SSVEP = pop_select( EEG_SSVEP,'notrial',t.tr2discard );
%figure; pop_plottopo(EEG_SSVEP, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);

% extract SSVEP
% for trials
t.fftrate =EEG.srate*5;

FFTData_tr=nan(size(EEG_SSVEP.data,3),numel(F.EEGChans),t.fftrate);
for i_ep = 1:size(EEG_SSVEP.data,3)
    FFTData_tr(i_ep,:,:)=abs(fft(EEG_SSVEP.data(1:64,:,i_ep),t.fftrate,2))*2/size(EEG_SSVEP.data,2);
end
FFTData_mean=abs(fft(mean(EEG_SSVEP.data(1:64,:,:),3),t.fftrate,2))*2/size(EEG_SSVEP.data,2);
xScale = ((0:size(FFTData_tr,3)-1)/size(FFTData_tr,3)) * EEG.srate;

% topoplot
[~, t.pos]=min(abs(xScale-14.16667));
figure;
colormap('hot')
topoplot( FFTData_mean(:,t.pos), EEG_SSVEP.chanlocs, 'chaninfo', EEG_SSVEP.chaninfo, ...
    'shading', 'flat', 'numcontour', 0,'conv','on', 'maplimits',[0 max(FFTData_mean(:,t.pos))],'plotchans', F.EEGChans);
title('SSVEP - averaged in time domain')
colorbar

figure;
colormap('hot')
topoplot( squeeze(mean(FFTData_tr(:,:,t.pos),1)), EEG_SSVEP.chanlocs, 'chaninfo', EEG_SSVEP.chaninfo, ...
    'shading', 'flat', 'numcontour', 0,'conv','on', 'maplimits',[0 max(squeeze(mean(FFTData_tr(:,:,t.pos),1)))],'plotchans', F.EEGChans);
title('SSVEP - averaged in frequency domain')
colorbar

% plot spectra
pl.elec2plot={'Oz'};
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({EEG_SSVEP.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
figure;
plot(xScale,FFTData_mean(pl.elec2plot_i,:),'r')
hold on
plot(xScale,squeeze(mean(FFTData_tr(:,pl.elec2plot_i,:),1)),'b')
xlabel('frequency in Hz')
ylabel('amplitude in \muV')
legend({'averaged in time domain','averaged in frequency domain'})
xlim([0 55])


%% TFA power around button press
% index block length
t.i1=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{1}); 
t.i2=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{2});
%.times = [cell2mat({EEG.event(t.i1).latency})' cell2mat({EEG.event(t.i2).latency})'];
t.eptime=307.17666;

EEG_exp = pop_epoch( EEG, F.ExpEpochTrigger(1), [0 t.eptime], 'epochinfo', 'yes');


t.bl=1;
t.frames=85;
t.triggerfreq = 6; % every 6 seconds
t.flickframes=[1 1 1 0 0 0];
for i_bl = 1:6
    % trigger times for
    t.t_trig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==15}});
    t.t_btrig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==30}});
    
    %exptiming.frametime{i_bl}=  linspace(t.t_trig(1),t.t_trig(end),((numel(t.t_trig)-1)*t.frames*t.triggerfreq)+1);
    t.frametimes=[];
    for i_tr = 1:numel(t.t_trig)-1
        t.frametimes=[ t.frametimes linspace(t.t_trig(i_tr),t.t_trig(i_tr+1),t.frames*t.triggerfreq+1)];
    end
    exptiming.frametime{i_bl}=sort(unique(t.frametimes));
    
    % create vector with on of cycles
    exptiming.framemat{i_bl}=   [repmat(t.flickframes,1,(numel(exptiming.frametime{i_bl})-1)/numel(t.flickframes)) 1];
    
    % create vector with trigger onset
    exptiming.triggermat{i_bl}=    zeros(1,numel(exptiming.framemat{i_bl}));
    exptiming.triggers{i_bl}=  t.t_trig(t.t_trig>=exptiming.frametime{i_bl}(1)&t.t_trig<=exptiming.frametime{i_bl}(end));
    exptiming.triggermat{i_bl}(arrayfun(@(x) find(exptiming.frametime{i_bl}<=x,1,'last'), exptiming.triggers{i_bl}))=1;
    
    % create vector with response onsets
    exptiming.respmat{i_bl}=    zeros(1,numel(exptiming.framemat{i_bl}));
    exptiming.responses{i_bl}=  t.t_btrig(t.t_btrig>=exptiming.frametime{i_bl}(1)&t.t_btrig<=exptiming.frametime{i_bl}(end));
    exptiming.respmat{i_bl}(arrayfun(@(x) find(exptiming.frametime{i_bl}<=x,1,'last'), exptiming.responses{i_bl}))=1;
    
    % function to get only first output: subsref(test(input),struct('type','()','subs',{{1}}))
%     t.shifts = arrayfun(@(x) subsref(findstr(exptiming.framemat{1,1}(x-numel(t.flickframes)+1:end),[1 1]),struct('type','()','subs',{{1}}))-numel(t.flickframes),...
%         find(exptiming.respmat{i_bl}==1));
    t.shifts = arrayfun(@(x) subsref(findstr(exptiming.framemat{i_bl}(x-round(numel(t.flickframes)/2):end),...
        t.flickframes(1:end/2)),struct('type','()','subs',{{1}}))-round(numel(t.flickframes)/2)-1,...
        find(exptiming.respmat{i_bl}==1));
    exptiming.respmat_sh{i_bl}=zeros(size(exptiming.respmat{i_bl}));
    exptiming.respmat_sh{i_bl}(find(exptiming.respmat{i_bl})+t.shifts)=1;
    
    % t.tt=[exptiming.frametime{i_bl}; exptiming.framemat{i_bl};  exptiming.triggermat{i_bl}; exptiming.respmat{i_bl}; exptiming.respmat_sh{i_bl}];
end

for i_bl = 1:6
    t.ev=find(exptiming.respmat_sh{i_bl});
    for i_ev = 1:numel(t.ev)
        EEG_exp = pop_editeventvals(EEG_exp,'append',{1 [] [] [] 1},...
            'changefield',{2 'epoch' i_bl},...
            'changefield',{2 'latency' exptiming.frametime{i_bl}(t.ev(i_ev))},...
            'changefield',{2 'type' 35});
    end
end
% pop_eegplot(EEG_exp,1,1,1)
EEG_expep = pop_epoch( EEG_exp, {'35'}, [-2.5 2.5], 'epochinfo', 'yes');
EEG_expep = pop_epoch( EEG_exp, {'30'}, [-2.5 2.5], 'epochinfo', 'yes');
EEG_expep = pop_rmbase( EEG_expep, []);
EEG_expep = eegF_Detrend(EEG_expep);

%filter
EEG_expep = pop_eegfiltnew(EEG_expep, 0.5, [], EEG_expep.pnts, 0, [], 0);

%pop_eegplot(EEG_expep,1,1,1)
%figure; pop_plottopo(EEG_expep, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);


% calculate induced TFA (average in frequency domain)
F.TFAfreqs=[5:0.5:30];
TFA.data=nan(numel(F.TFAfreqs), EEG_expep.pnts, numel(F.EEGChans),2);
fprintf(1,'calculating %1.0f TFAs: ',numel(F.EEGChans))
for i_el=1:numel(F.EEGChans) %loop electrodes
    fprintf(1,'%1.0f',mod(F.EEGChans(i_el),10))
    [TFA.data(:,:,i_el,1),TFA.t,TFA.f]=traces2TFA_4(squeeze(EEG_expep.data(i_el,:,:)),F.TFAfreqs,EEG_expep.srate,5);
    [TFA.data(:,:,i_el,2),TFA.t,TFA.f]=traces2TFA_4(squeeze(nanmean(EEG_expep.data(i_el,:,:),3))',F.TFAfreqs,EEG_expep.srate,5);
end
fprintf(1,'...done\n')
TFA.t=EEG_expep.times;

% baseline correct data
TFA.baseline=[-2250 -1750];
TFA.databasecorr = bsxfun(@minus, TFA.data, mean(TFA.data(:,eeg_time2points(TFA.baseline(1),TFA.t):eeg_time2points(TFA.baseline(2),TFA.t),:,:),2));

%% plot data
% pl.elec2plot = {'Oz'};
% pl.elec2plot = {'PO3'};
% pl.elec2plot = {'PO4'};
% pl.elec2plot = {'PO7'};
pl.elec2plot = {'PO8'};
% pl.elec2plot = {'C3'};
% pl.elec2plot = {'Cz'};
pl.xlims=[-2000 2000];
t.elec2plot=find(strcmp({EEG_expep.chanlocs.labels},pl.elec2plot));
[t.t t.ind1]=min(abs(TFA.t-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.t-pl.xlims(2)));

% raw
t.clims1=[0 1]*max(max(max(TFA.data(:,t.ind1:t.ind2,t.elec2plot,1))));
t.clims2=[0 1]*max(max(max(TFA.data(:,t.ind1:t.ind2,t.elec2plot,2))));
figure;
subplot(2,1,1)
imagesc(TFA.t,TFA.f,TFA.data(:,:,t.elec2plot,1),t.clims1)
set(gca,'YDir','normal')
title(['tfa, for channel ' pl.elec2plot{1} '; averaged in frequency domain'])
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)

subplot(2,1,2)
imagesc(TFA.t,TFA.f,TFA.data(:,:,t.elec2plot,2),t.clims2)
set(gca,'YDir','normal')
title(['tfa, for channel ' pl.elec2plot{1} '; averaged in time domain'])
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)

% basecorr
t.clims1=[-1 1]*max(max(max(abs(TFA.databasecorr(:,t.ind1:t.ind2,t.elec2plot,1)))));
t.clims2=[-1 1]*max(max(max(abs(TFA.databasecorr(:,t.ind1:t.ind2,t.elec2plot,2)))));
figure;
subplot(2,1,1)
imagesc(TFA.t,TFA.f,TFA.databasecorr(:,:,t.elec2plot,1),t.clims1)
set(gca,'YDir','normal')
title(['corrected tfa, for channel ' pl.elec2plot{1} '; averaged in frequency domain'])
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)

subplot(2,1,2)
imagesc(TFA.t,TFA.f,TFA.databasecorr(:,:,t.elec2plot,2),t.clims2)
set(gca,'YDir','normal')
title(['corrected tfa, for channel ' pl.elec2plot{1} '; averaged in time domain'])
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)


% plot timecourse of frequency
pl.freq2plot=14.16667;
% pl.freq2plot=14.16667*2;
[t.t t.ind3]=min(abs(TFA.f-pl.freq2plot));

figure;
subplot(2,1,1)
plot(TFA.t,TFA.databasecorr(t.ind3,:,t.elec2plot,1))
ylim([-1.2 1.2]*max(abs(TFA.databasecorr(t.ind3,t.ind1:t.ind2,t.elec2plot,1))))
xlabel('time in ms')
ylabel('amplitude in \muV')
xlim(pl.xlims)
title(sprintf('timecourse signal %1.2fHz; averaged in time domain\n channel %s; baseline [%1.0f %1.0f]ms',...
    pl.freq2plot, pl.elec2plot{1},TFA.baseline))
grid on

subplot(2,1,2)
plot(TFA.t,TFA.databasecorr(t.ind3,:,t.elec2plot,2))
ylim([-1.2 1.2]*max(abs(TFA.databasecorr(t.ind3,t.ind1:t.ind2,t.elec2plot,2))))
xlabel('time in ms')
ylabel('amplitude in \muV')
xlim(pl.xlims)
title(sprintf('timecourse signal %1.2fHz; averaged in frequency domain\n channel %s; baseline [%1.0f %1.0f]ms',...
    pl.freq2plot, pl.elec2plot{1},TFA.baseline))
grid on











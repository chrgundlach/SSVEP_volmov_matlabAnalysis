%% script in order to calculate time frequency analysis of main experiment



clear all

%% parameters
F.PathIn                = 'F:\work\data\SSVEP_volmov\EEG\RAW';
F.PathOut               = 'F:\work\data\SSVEP_volmov\EEG\TFA_SCADS';
F.Subjects2Use          = [1:5];
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.VEOGChans             = [];
F.HEOGChans             = [];
F.ChanlocFile           = 'C:\Users\HP-User\matlab\Skripte\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.RespTrigger           = {'30'};
F.RespEpoch             = [-3.5 3.5]; % in s
F.SSVEPTrigger          = {'15'};

F.ExpEpochTrigger       = {'11','12'};
F.eptime                = 359.5; % time of each block in s
F.Blocks                = 6; % number of blocks

F.FrameRate             = 85;
F.TriggerRate           = 6; % SSVEP trigger frequency (e.g. 6 = every 6s)
F.flickframes           = [1 1 1 0 0 0]; % on-off frames for SSVEP

F.TFAfreqs              = [4:0.2:45];
    

%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp.bdf ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    EEG=pop_biosig(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)));
    % pop_eegplot(EEG,1,1,1);
    
    % chanlocs
    EEG.chanlocs=pop_chanedit(EEG.chanlocs,'load',{F.ChanlocFile,'filetype','besa (elp)'}); % Load Channel Locations
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    
    % rereference
    EEG = pop_reref( EEG, [],'exclude',[F.EEGChans(end)+1:EEG.nbchan] ); % average
    % pop_eegplot(EEG,1,1,1,1)
    
    %% TFA power around button press
    % induced and evoked...requires shift of button press to SSVEP cycle beginning
    
    % index block length
    t.i1=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{1});
    t.i2=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{2});
    %.times = [cell2mat({EEG.event(t.i1).latency})' cell2mat({EEG.event(t.i2).latency})'];
    
    EEG_exp = pop_epoch(EEG, F.ExpEpochTrigger(1), [0 F.eptime], 'epochinfo', 'yes');
    % pop_eegplot(EEG_exp,1,1,1,1)
    
    %% loop for blocks
    % create time vectors for visual stimulation
    t.bl=1;
    for i_bl = 1:F.Blocks
        % trigger times for
        t.t_trig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.SSVEPTrigger{1})}}); % SSVEP trigger
        t.t_btrig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.RespTrigger{1})}}); % behavior trigger
        
        t.frametimes=[];
        % create time for each frame (fit between triggers)
        for i_tr = 1:numel(t.t_trig)-1
            t.frametimes=[ t.frametimes linspace(t.t_trig(i_tr),t.t_trig(i_tr+1),F.FrameRate*F.TriggerRate+1)];
        end
        % comepensate for missing last SSVEP trigger...add time towards the end
        t.frametimes=[ t.frametimes (0:median(diff(t.frametimes)):(median(diff(t.t_trig))))+t.frametimes(end)];
        t.frametimes=sort(unique(t.frametimes));
        % discard frames outside epoch
        exptiming.frametime{i_bl}=t.frametimes(t.frametimes<EEG_exp.times(end));
        
        % create vector with on-off cycles
        exptiming.framemat{i_bl}=[repmat(F.flickframes,1,floor((numel(exptiming.frametime{i_bl}))/numel(F.flickframes)))];
        exptiming.framemat{i_bl}=[exptiming.framemat{i_bl} F.flickframes(1:(numel(exptiming.frametime{i_bl})-numel(exptiming.framemat{i_bl})))];
        
        % create vector with trigger onset
        exptiming.triggermat{i_bl}= zeros(1,numel(exptiming.framemat{i_bl}));
        exptiming.triggers{i_bl}=  t.t_trig(t.t_trig>=exptiming.frametime{i_bl}(1)&t.t_trig<=exptiming.frametime{i_bl}(end));
        exptiming.triggermat{i_bl}(arrayfun(@(x) find(exptiming.frametime{i_bl}<=x,1,'last'), exptiming.triggers{i_bl}))=1;
        
        % create vector with response onsets
        exptiming.respmat{i_bl}=    zeros(1,numel(exptiming.framemat{i_bl}));
        exptiming.responses{i_bl}=  t.t_btrig(t.t_btrig>=exptiming.frametime{i_bl}(1)&t.t_btrig<=exptiming.frametime{i_bl}(end));
        exptiming.respmat{i_bl}(arrayfun(@(x) find(exptiming.frametime{i_bl}<=x,1,'last'), exptiming.responses{i_bl}))=1;
        
        % function to get only first output: subsref(test(input),struct('type','()','subs',{{1}}))
        %     t.shifts = arrayfun(@(x) subsref(findstr(exptiming.framemat{1,1}(x-numel(F.flickframes)+1:end),[1 1]),struct('type','()','subs',{{1}}))-numel(F.flickframes),...
        %         find(exptiming.respmat{i_bl}==1));
        t.shifts = arrayfun(@(x) subsref(findstr(exptiming.framemat{i_bl}(x-round(numel(F.flickframes)/2):end),...
            F.flickframes(1:end/2)),struct('type','()','subs',{{1}}))-round(numel(F.flickframes)/2)-1,...
            find(exptiming.respmat{i_bl}==1));
        exptiming.respmat_sh{i_bl}=zeros(size(exptiming.respmat{i_bl}));
        exptiming.respmat_sh{i_bl}(find(exptiming.respmat{i_bl})+t.shifts)=1;
        
        % t.tt=[exptiming.frametime{i_bl}; exptiming.framemat{i_bl};  exptiming.triggermat{i_bl}; exptiming.respmat{i_bl}; exptiming.respmat_sh{i_bl}];
    end
    
    %% write new events
    for i_bl = 1:6
        t.ev=find(exptiming.respmat_sh{i_bl});
        for i_ev = 1:numel(t.ev)
            EEG_exp = pop_editeventvals(EEG_exp,'append',{1 [] [] [] 1},...
                'changefield',{2 'epoch' i_bl},...
                'changefield',{2 'latency' exptiming.frametime{i_bl}(t.ev(i_ev))},...
                'changefield',{2 'type' 35});
        end
    end
    % resample
    EEG_exp = pop_resample(EEG_exp, 250);
    %% extract new epocs
    % pop_eegplot(EEG_exp,1,1,1)
    EEG_expep_sh = pop_epoch( EEG_exp, {'35'}, F.RespEpoch, 'epochinfo', 'yes');
    EEG_expep = pop_epoch( EEG_exp, {'30'}, F.RespEpoch, 'epochinfo', 'yes');
    %EEG_expep = pop_epoch( EEG_exp, {'15'}, F.RespEpoch, 'epochinfo', 'yes');
    
    EEG_expep = eegF_Detrend(EEG_expep);
    EEG_expep_sh = eegF_Detrend(EEG_expep_sh);
    
    %filter
    TFA.params.filter = {[] []};
    %TFA.params.filter = {0.5 []};
    %EEG_expep = pop_eegfiltnew(EEG_expep, TFA.params.filter{1}, TFA.params.filter{2}, EEG_expep.pnts, 0, [], 0);
    %EEG_expep_sh = pop_eegfiltnew(EEG_expep_sh, TFA.params.filter{1}, TFA.params.filter{2}, EEG_expep_sh.pnts, 0, [], 0);
    
    
    %pop_eegplot(EEG_expep,1,1,1)
    %figure; pop_plottopo(EEG_expep, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    %figure; pop_plottopo(EEG_expep_sh, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    
    TFA.electrodes = EEG_expep.chanlocs;
    TFA.params.srate = EEG_exp.srate;
    
    %TFA.wave_cycles=7;
    TFA.wave_cycles=[4 60];
    TFA.wave_cycles=[3.5 (F.TFAfreqs(end))/(F.TFAfreqs(1)/3.5)];
    
    %% calculate induced TFA (average in frequency domain) and evoked TFA
    TFA.data_induced=nan(numel(F.TFAfreqs), EEG_expep.pnts/2, numel(F.EEGChans));
    TFA.data_evoked=TFA.data_induced;
    TFA.alltrials = EEG_expep.trials;
    
    fprintf(1,'calculating %1.0f TFAs: ',numel(F.EEGChans))
    for i_el=1:numel(F.EEGChans) %loop electrodes
        fprintf(1,'%1.0f',mod(F.EEGChans(i_el),10))
        % preliminary artifact correction
        t.data = squeeze(EEG_expep.data(i_el,:,:));
        TFA.trials_induced{i_el}=find(sum((t.data>(mean(t.data(:))+6*std(t.data(:))))|(t.data<(mean(t.data(:))-6*std(t.data(:)))))==0);
        t.data = squeeze(EEG_expep_sh.data(i_el,:,:));
        TFA.trials_evoked{i_el}=find(sum((t.data>(mean(t.data(:))+6*std(t.data(:))))|(t.data<(mean(t.data(:))-6*std(t.data(:)))))==0);
        
        
        %[t.tfadata,TFA.t,TFA.f] = traces2TFA_4(squeeze(EEG_expep.data(i_el,:,TFA.trials_induced{i_el})),F.TFAfreqs,EEG_expep.srate,TFA.wave_cycles);
        [t.tfadata,TFA.t,TFA.f,t.D] = mwavelTFA(squeeze(EEG_expep.data(i_el,:,TFA.trials_induced{i_el})),F.TFAfreqs,EEG_expep.srate,TFA.wave_cycles,1,3.5);
        TFA.data_induced(:,:,i_el) = imresize(mean(t.tfadata,3),[numel(F.TFAfreqs),EEG_expep.pnts/2]);
        
        %[t.tfadata,TFA.t,TFA.f] = traces2TFA_4(squeeze(nanmean(EEG_expep_sh.data(i_el,:,TFA.trials_evoked{i_el}),3))',F.TFAfreqs,EEG_expep.srate,TFA.wave_cycles);
        [t.tfadata,TFA.t,TFA.f,TFA.D] = mwavelTFA(squeeze(nanmean(EEG_expep_sh.data(i_el,:,TFA.trials_evoked{i_el}),3))',F.TFAfreqs,EEG_expep.srate,TFA.wave_cycles,1,3.5);
        TFA.data_evoked(:,:,i_el)=imresize(t.tfadata,[numel(F.TFAfreqs),EEG_expep_sh.pnts/2]);
        % figure; imagesc(TFA.t,TFA.f,mean(t.tfadata,3)); set(gca,'ydir','normal')
    end
    fprintf(1,'...done\n')
    TFA.t=imresize(EEG_expep.times,[1,EEG_expep_sh.pnts/2]);
    
    %% save
    TFA.savetime=datestr(now);
    fprintf(1,'|| saving file ||  %s\\VP%02.0f_exp_tfa.mat ||\n', F.PathOut,F.Subjects2Use(i_sub)')
    save(sprintf('%s\\VP%02.0f_exp_tfa.mat', F.PathOut,F.Subjects2Use(i_sub)), 'TFA')   
    
end
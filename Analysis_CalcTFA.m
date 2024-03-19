%% script in order to calculate time frequency analysis of main experiment



clear all

%% parameters
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\RAW';
F.PathInEOG             = 'D:\work\data\SSVEP_volmov\EEG\EOG_Data\';
%F.PathOut               = 'F:\work\data\SSVEP_volmov\EEG\TFA_SCADS_Gabor';
F.PathOut               = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_noblinks';
F.Subjects2Use          = [1:20];
%F.Subjects2Use          = [1];
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.VEOGChans             = [];
F.HEOGChans             = [];
F.ChanlocFile           = 'C:\Users\psy05cvd\Dropbox\work\matlab\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.RespTrigger           = {'30'};
F.RespEpoch             = [-5.5 5.5]; % in s
F.SSVEPTrigger          = {'15'};

F.ExpEpochTrigger       = {'11','12'};
F.eptime                = 359.5; % time of each block in s
F.Blocks                = 6; % number of blocks

F.FrameRate             = 85;
F.TriggerRate           = 6; % SSVEP trigger frequency (e.g. 6 = every 6s)
F.flickframes           = [1 1 1 0 0 0]; % on-off frames for SSVEP

F.TFAfreqs              = [4-(1/3):0.25:45];
F.TFAFlag               = [2]; % 1=wavelet; 2=gabor    

%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp.bdf ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
%     EEG=pop_biosig(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)));
    EEG=pop_readbdf(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)),[],73,[]);
    % pop_eegplot(EEG,1,1,1);
    
    % load EOG data
    EOG = open(sprintf('%s\\VP%02.0f_exp_EOG.mat',F.PathInEOG,F.Subjects2Use(i_sub)));
    
    % rereference
    EEG = pop_reref( EEG, [],'exclude',[F.EEGChans(end)+1:EEG.nbchan] ); % average    
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
    EEG_exp=eeg_checkset(EEG_exp);
    
    %% extract new epocs + artifact correction
    % pop_eegplot(EEG_exp,1,1,1)
    
    % adjusted for artifact correction + SCADS implementation
    % resample
    EEG_exp = pop_resample(EEG_exp, 256);
    
    % chanlocs
    EEG_exp.chanlocs=pop_chanedit(EEG_exp.chanlocs,'load',{F.ChanlocFile,'filetype','besa (elp)'}); % Load Channel Locations
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    EEG_ext = pop_select(EEG_exp,'channel',65:72);
    EEG_exp = pop_select(EEG_exp,'channel',1:64);
    
    % epoch (depending on F.TFAFlag)
    switch F.TFAFlag
        case 1
            t.ep_time = F.RespEpoch;
        case 2
            t.ep_time = F.RespEpoch+[-0.5 0.5];
    end
    EEG_expep_n = pop_epoch( EEG_exp, {'30'}, t.ep_time, 'epochinfo', 'yes');
    EEG_expep_sh = pop_epoch( EEG_exp, {'35'}, t.ep_time, 'epochinfo', 'yes');
    EEG_ext = pop_epoch( EEG_ext, {'30'}, t.ep_time, 'epochinfo', 'yes');
    
    % get timing
    % get timing (adjusted)
    art.resp_time_n=[];
    for i_tr = 1:EEG_expep_n.trials
        t.ix=EEG_expep_n.epoch(i_tr).eventurevent{dsearchn(cell2mat(EEG_expep_n.epoch(i_tr).eventlatency)',0)};
        art.resp_time_n(i_tr)=EEG_expep_n.urevent(t.ix).latency./EEG_expep_n.srate;
    end
    art.resp_time_sh=[];
    for i_tr = 1:EEG_expep_sh.trials
        t.ix=EEG_expep_sh.epoch(i_tr).eventurevent{dsearchn(cell2mat(EEG_expep_sh.epoch(i_tr).eventlatency)',0)};
        art.resp_time_sh(i_tr)=EEG_expep_sh.urevent(t.ix).latency./EEG_expep_sh.srate;
    end
    
    % discard trials with blinks around movement period [-2 +2] s
    t.rejindex_n = false(1,EEG_expep_n.trials);
    t.rejindex_sh = false(1,EEG_expep_sh.trials);
    % index trials in certain time interval
    art.win_blinks = [-3000 3000];
    for i_tr = 1:EEG_expep_n.trials
        if EOG.EOG.art.blinknum_in_resptrials_n(i_tr)
            t.rejindex_n(i_tr)=any((art.win_blinks(1) <= EOG.EOG.art.blink_time_in_resptrials_n{i_tr}) &...
                (art.win_blinks(2) >= EOG.EOG.art.blink_time_in_resptrials_n{i_tr}));
        end
    end
    for i_tr = 1:EEG_expep_sh.trials
        if EOG.EOG.art.blinknum_in_resptrials_sh(i_tr)
            t.rejindex_sh(i_tr)=any((art.win_blinks(1) <= EOG.EOG.art.blink_time_in_resptrials_sh{i_tr}) &...
                (art.win_blinks(2) >= EOG.EOG.art.blink_time_in_resptrials_sh{i_tr}));
        end
    end
    % some later errors with eeglab if trialnumber < 2
    if diff([sum(t.rejindex_n) numel(t.rejindex_n)])<2
        t.rejindex_n([1 2])=false;
    end
    if diff([sum(t.rejindex_sh) numel(t.rejindex_sh)])<2
        t.rejindex_sh([1 2])=false;
    end
    
    
    % detrend
    EEG_expep_n = eegF_Detrend(EEG_expep_n,t.ep_time);
    EEG_expep_sh = eegF_Detrend(EEG_expep_sh,t.ep_time);
    EEG_ext = eegF_Detrend(EEG_ext,t.ep_time);
    %pop_eegplot(EEG_expep_n,1,1,1)
    %pop_eegplot(EEG_expep_sh,1,1,1)
    %figure; plot(EEG_ext.times,mean(EEG_ext.data,3)); legend({EEG_ext.chanlocs.labels})
    %figure; plot(EEG_ext.times,std(EEG_ext.data,1,3)); legend({EEG_ext.chanlocs.labels})
    
    % run SCADS to statistically detetct artifacts for unshifted trials
    %[EEG_expep_n art.SCADS_Trials2Del_n art.SCADS_SumData_n] = SCADS_Pass1(EEG_expep_n,F.RespEpoch,[6 6 15],1);
    % discard trials indexed as containing artifacts
%     t.rejindex = EEG_expep_n.reject.rejmanual;
%     EEG_expep_n = pop_rejepoch(EEG_expep_n, EEG_expep_n.reject.rejmanual, 0);
    
    % run SCADS to statistically detetct artifacts for shifted trials
    %[EEG_expep_sh art.SCADS_Trials2Del_sh art.SCADS_SumData_sh] = SCADS_Pass1(EEG_expep_sh,F.RespEpoch,[6 6 15],1);
    % discard trials indexed as containing artifacts
%     EEG_expep_sh = pop_rejepoch(EEG_expep_sh, EEG_expep_sh.reject.rejmanual, 0);
    
    % reject trials for external signals
    EEG_ext_rej = EEG_ext;
    %EEG_ext_rej = pop_rejepoch(EEG_ext, t.rejindex, 0);
    %figure; plot(EEG_ext_rej.times,mean(EEG_ext_rej.data,3)); legend({EEG_ext_rej.chanlocs.labels})
    %figure; plot(EEG_ext_rej.times,std(EEG_ext_rej.data,1,3)); legend({EEG_ext_rej.chanlocs.labels})
   
    % rereference
    %EEG_expep_n=pop_reref(EEG_expep_n,[],'refstate',0);
    %EEG_expep_sh=pop_reref(EEG_expep_sh,[],'refstate',0);
    %pop_eegplot(EEG_expep_n,1,1,1)
    %pop_eegplot(EEG_expep_sh,1,1,1)
    
   
    %pop_eegplot(EEG_expep,1,1,1)
    %figure; pop_plottopo(EEG_expep, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    %figure; pop_plottopo(EEG_expep_sh, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    
    
    %% reject trials with artifacts
    % reject trials
    EEG_expep_n = pop_rejepoch(EEG_expep_n, t.rejindex_n, 0);
    EEG_expep_sh = pop_rejepoch(EEG_expep_sh, t.rejindex_sh, 0);
    EEG_ext = pop_rejepoch(EEG_ext, t.rejindex_n, 0);
    
    EEG_surrogate =  pop_select( EEG_expep_n,'time',F.RespEpoch);
        
    %% calculate induced TFA (average in frequency domain) and evoked TFA
    TFA.electrodes = EEG_expep_n.chanlocs;
    TFA.params.srate = EEG_expep_n.srate;
    
    %TFA.wave_cycles=7;
    TFA.wave_cycles=[4 60];
    TFA.wave_cycles=[4 (F.TFAfreqs(end))/(F.TFAfreqs(1)/4)];
    TFA.wave_cycles=[6 (F.TFAfreqs(end))/(F.TFAfreqs(1)/6)];
    
    %filter
    TFA.params.filter = {[] []};
    %TFA.params.filter = {0.5 []};
    %EEG_expep_n = pop_eegfiltnew(EEG_expep_n, TFA.params.filter{1}, TFA.params.filter{2}, EEG_expep_n.pnts, 0, [], 0);
    %EEG_expep_sh = pop_eegfiltnew(EEG_expep_sh, TFA.params.filter{1}, TFA.params.filter{2}, EEG_expep_sh.pnts, 0, [], 0);
    
    TFA.data_induced=nan(numel(F.TFAfreqs), EEG_surrogate.pnts/2, numel(F.EEGChans));
    TFA.data_evoked=TFA.data_induced;
    TFA.alltrials = EEG_expep_n.trials;
    
    
    switch F.TFAFlag
        case 1
            fprintf(1,'calculating %1.0f TFAs: ',numel(F.EEGChans))
            TFA.type = 'wavelet';
            for i_el=1:numel(F.EEGChans) %loop electrodes
                fprintf(1,'%1.0f',mod(F.EEGChans(i_el),10))
                TFA.type = 'wavelet';
                
                [t.tfadata,TFA.t,TFA.f,t.D] = mwavelTFA(squeeze(EEG_expep_n.data(i_el,:,:)),F.TFAfreqs,EEG_expep_n.srate,TFA.wave_cycles,1,3.5);
                TFA.data_induced(:,:,i_el) = imresize(mean(t.tfadata,3),[numel(F.TFAfreqs),EEG_expep_n.pnts/2]);
                
                [t.tfadata,TFA.t,TFA.f,TFA.D] = mwavelTFA(squeeze(nanmean(EEG_expep_sh.data(i_el,:,:),3))',F.TFAfreqs,EEG_expep_sh.srate,TFA.wave_cycles,1,3.5);
                TFA.data_evoked(:,:,i_el)=imresize(t.tfadata,[numel(F.TFAfreqs),EEG_expep_sh.pnts/2]);
                % figure; imagesc(TFA.t,TFA.f,mean(t.tfadata,3)); set(gca,'ydir','normal')
            end
        case 2
            TFA.f = F.TFAfreqs;
            fprintf(1,'calculating %1.0f TFAs: \n',numel(TFA.f))
            TFA.type = 'gabor';
            TFA.params.gabor_FWHM_freq = 0.5;
            for i_freq = 1:numel(TFA.f)
                % induced
                EEG_Gabor_gab = eegF_Gabor(EEG_expep_n, TFA.f(i_freq),TFA.params.gabor_FWHM_freq);
                EEG_Gabor_gab = pop_select( EEG_Gabor_gab,'time',F.RespEpoch);
                TFA.data_induced(i_freq,:,:)=imresize(mean(EEG_Gabor_gab.data,3),[EEG_Gabor_gab.nbchan,EEG_Gabor_gab.pnts/2])';
                
                % evoked
                EEG_t = EEG_expep_sh; EEG_t=pop_select(EEG_t,'trial',1); EEG_t.data=mean(EEG_expep_sh.data,3);
                EEG_Gabor_gab = eegF_Gabor(EEG_t, TFA.f(i_freq),TFA.params.gabor_FWHM_freq);
                EEG_Gabor_gab = pop_select( EEG_Gabor_gab,'time',F.RespEpoch);
                TFA.data_evoked(i_freq,:,:)=imresize(mean(EEG_Gabor_gab.data,3),[EEG_Gabor_gab.nbchan,EEG_Gabor_gab.pnts/2])';
            end
    end
    fprintf(1,'...done\n')
    TFA.t=imresize(EEG_surrogate.times,[1,EEG_surrogate.pnts/2]);
    
    %% save
    TFA.art=art;
    TFA.savetime=datestr(now);
    fprintf(1,'|| saving file ||  %s\\VP%02.0f_exp_tfa.mat ||\n', F.PathOut,F.Subjects2Use(i_sub)')
    save(sprintf('%s\\VP%02.0f_exp_tfa.mat', F.PathOut,F.Subjects2Use(i_sub)), 'TFA')   
    
end
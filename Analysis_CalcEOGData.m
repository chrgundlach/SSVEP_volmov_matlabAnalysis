%% script to analyse EOG data
% in relation to button presses
% look at SSVEPs around blinks
% different RESS settings for VP 17 


clearvars

%% parameters
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\RAW';
F.PathOut               = 'D:\work\data\SSVEP_volmov\EEG\EOG_Data';
% F.Subjects2Use          = [1:20];
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27];
% F.Subjects2Use          = [19];
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.VEOGChans             = [65 66];
F.HEOGChans             = [67 68];
F.ChanlocFile           = 'C:\Users\psy05cvd\Dropbox\work\matlab\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.RespTrigger           = {'30'};
F.RespEpoch             = [-6 6]; % in s
F.EOGEpoch              = [-5.5 5.5]; % in s
F.SSVEPTrigger          = {'15'};
F.BlinkTrigger          = {'50'};
F.CSD_flag              = 1;

F.ExpEpochTrigger       = {'11','12'};
F.eptime                = 360; % time of each block in s
F.Blocks                = 6; % number of blocks

F.FrameRate             = 85;
F.TriggerRate           = 6; % SSVEP trigger frequency (e.g. 6 = every 6s)
F.flickframes           = [1 1 1 0 0 0]; % on-off frames for SSVEP

F.TFAfreqs              = [2-(1/3):0.25:31];
F.TFAFlag               = [2]; % 1=wavelet; 2=gabor
F.CSD_flag              = 1; % 0 = no; 1 = yes




%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp.bdf ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
%     EEG=pop_biosig(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)));
    
    if F.Subjects2Use(i_sub) ~= 24
        EEG=pop_readbdf(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)),[],F.EMGChans(end)+1,[]);
    else
        EEG = pop_loadset('filename','VP24_exp.set','filepath',F.PathIn);
    end
    % troubleshooting for subject 24
%     EEG = pop_loadset('filename','VP24_exp.set','filepath','F:\\work\\data\\SSVEP_volmov\\EEG\\RAW\\');
%     t.idx = strcmp({EEG.event.type},'boundary');
%     EEG = pop_editeventvals(EEG,'delete',find(t.idx));
%     for i_ev = 1:numel(EEG.event)
%         EEG.event(i_ev).type = str2double(EEG.event(i_ev).type);
%     end
    
    % pop_eegplot(EEG,1,1,1);
    %numel(find(cell2mat({EEG.event.type})==15))
    
    % rereference
    EEG = pop_reref( EEG, [],'exclude',[F.EEGChans(end)+1:EEG.nbchan] ); % average
    % pop_eegplot(EEG,1,1,1,1)
    
    EEG = pop_resample(EEG, 256);
    %% epoch around button presses
    % requires shift of button press to SSVEP cycle beginning
    
    % index block length
%     t.i1=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{1});
%     t.i2=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{2});
    %.times = [cell2mat({EEG.event(t.i1).latency})' cell2mat({EEG.event(t.i2).latency})'];
    
    EEG_exp = pop_epoch(EEG, F.ExpEpochTrigger(1), [0 F.eptime], 'epochinfo', 'yes');
    if EEG_exp.trials < 6
        EEG_exp = pop_epoch(EEG, F.ExpEpochTrigger(1), [0 F.eptime-0.5], 'epochinfo', 'yes');
    end
    % pop_eegplot(EEG_exp,1,1,1,1)
    
    % do some checking
    %     checking.trials(i_sub) = EEG_exp.trials;
    %     checking.timing_first{i_sub}=...
    %         cell2mat(cellfun(@(x) EEG_exp.epoch(x).eventlatency(find(cell2mat(EEG_exp.epoch(x).eventtype)==15,1,'first')), num2cell(1:EEG_exp.trials)));
    %     checking.timing_last{i_sub}=...
    %         cell2mat(cellfun(@(x) EEG_exp.epoch(x).eventlatency(find(cell2mat(EEG_exp.epoch(x).eventtype)==15,1,'last')), num2cell(1:EEG_exp.trials)));
    %     checking.timing_num{i_sub}=cellfun(@(x) numel(find(cell2mat(EEG_exp.epoch(x).eventtype)==15)), num2cell(1:EEG_exp.trials));
    %     cellfun(@(x) numel(find(cell2mat(EEG_exp.epoch(x).eventtype)==15)), num2cell(1:6))
    %     cellfun(@(x) EEG_exp.epoch(x).eventlatency(find(cell2mat(EEG_exp.epoch(x).eventtype)==15,1,'first')), num2cell(1:6),'UniformOutput',false)
    
    %% find blinks
    % bipolarize data
    EEG_exp_EOG=eegF_Bipolarize(EEG_exp);
    EEG_exp_EOG=pop_select(EEG_exp_EOG,'channel',{'VEOG','HEOG'});
    EEG_exp_EOG_f = pop_eegfiltnew(EEG_exp_EOG,0.5, 0, 8*EEG_exp_EOG.srate, 0, [], 0);
    %figure; plot(EEG_exp_EOG_f.times,EEG_exp_EOG_f.data(1,:,1));
    for i_bl = 1:F.Blocks
        [blinks.peaks blinks.peaklocs] = findpeaks(EEG_exp_EOG_f.data(1,:,i_bl),'MinPeakDistance',EEG_exp_EOG_f.srate/2,'MinPeakHeight',100);
        %figure; plot(EEG_exp_EOG_f.times,EEG_exp_EOG_f.data(1,:,i_bl)); 
        %hold on; plot(EEG_exp_EOG_f.times(blinks.peaklocs),EEG_exp_EOG_f.data(1,blinks.peaklocs,i_bl),'or')
        
        % write events of blinks
        for i_ev = 1:numel(blinks.peaklocs)
            EEG_exp = pop_editeventvals(EEG_exp,'append',{1 [] [] [] 1},...
                'changefield',{2 'epoch' i_bl},...
                'changefield',{2 'latency' blinks.peaklocs(i_ev)/EEG_exp_EOG_f.srate*1000},... % in seconds
                'changefield',{2 'type' F.BlinkTrigger{1}});
        end
    end
    EEG_exp=eeg_checkset(EEG_exp);
    % pop_eegplot(EEG_exp,1,1,1)
    
    %% loop for blocks
    % create time vectors for visual stimulation
    t.bl=1;
    for i_bl = 1:F.Blocks
        % trigger times for
        t.t_trig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.SSVEPTrigger{1})}}); % SSVEP trigger
        t.t_btrig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.RespTrigger{1})}}); % behavior trigger
        t.t_bltrig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.BlinkTrigger{1})}}); % blink trigger
        if isempty(t.t_trig)
            t.t_trig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{strcmp(EEG_exp.epoch(i_bl).eventtype,F.SSVEPTrigger{1})}}); % SSVEP trigger
        end
        if isempty(t.t_btrig)
            t.t_btrig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{strcmp(EEG_exp.epoch(i_bl).eventtype,F.RespTrigger{1})}}); % SSVEP trigger
        end
        if isempty(t.t_bltrig)
            t.t_bltrig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{strcmp(EEG_exp.epoch(i_bl).eventtype,F.BlinkTrigger{1})}}); % SSVEP trigger
        end
            
        
        t.frametimes=[];
        % create time for each frame (fit between triggers)
        for i_tr = 1:numel(t.t_trig)-1 % 1:numel(t.t_trig)-1
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
        exptiming.framemat{i_bl}(end+1)=1;
        
        % create vector with trigger onset
        exptiming.triggermat{i_bl}= zeros(1,numel(exptiming.framemat{i_bl}));
        exptiming.triggers{i_bl}=  t.t_trig(t.t_trig>=exptiming.frametime{i_bl}(1)&t.t_trig<=exptiming.frametime{i_bl}(end));
        exptiming.triggermat{i_bl}(arrayfun(@(x) find(exptiming.frametime{i_bl}<=x,1,'last'), exptiming.triggers{i_bl}))=1;
        
        % create vector with response onsets
        exptiming.respmat{i_bl}=    zeros(1,numel(exptiming.framemat{i_bl}));
        exptiming.responses{i_bl}=  t.t_btrig(t.t_btrig>=exptiming.frametime{i_bl}(1)&t.t_btrig<=exptiming.frametime{i_bl}(end));
        exptiming.respmat{i_bl}(arrayfun(@(x) find(exptiming.frametime{i_bl}<=x,1,'last'), exptiming.responses{i_bl}))=1;
        
        % create vector with blink onsets
        exptiming.blinkmat{i_bl}=    zeros(1,numel(exptiming.framemat{i_bl}));
        exptiming.blinks{i_bl}=  t.t_bltrig(t.t_bltrig>=exptiming.frametime{i_bl}(1)&t.t_bltrig<=exptiming.frametime{i_bl}(end));
        exptiming.blinkmat{i_bl}(arrayfun(@(x) find(exptiming.frametime{i_bl}<=x,1,'last'), exptiming.blinks{i_bl}))=1;
        
        % function to get only first output: subsref(test(input),struct('type','()','subs',{{1}}))
        %     t.shifts = arrayfun(@(x) subsref(findstr(exptiming.framemat{1,1}(x-numel(F.flickframes)+1:end),[1 1]),struct('type','()','subs',{{1}}))-numel(F.flickframes),...
        %         find(exptiming.respmat{i_bl}==1));
        % loopwise programming:
        % t.shifts2=[];
        % for i_resp = find(exptiming.respmat{i_bl}==1)
        %     t.idx = strfind(exptiming.framemat{i_bl}(i_resp-round(numel(F.flickframes)/2):end),F.flickframes(1:end/2));
        %     t.shifts2(end+1)=t.idx(1)-round(numel(F.flickframes)/2)-1;
        % end
        t.shifts = arrayfun(@(x) subsref(strfind(exptiming.framemat{i_bl}(x-round(numel(F.flickframes)/2):end),F.flickframes(1:end/2)),...
            struct('type','()','subs',{{1}}))-round(numel(F.flickframes)/2)-1,...
            find(exptiming.respmat{i_bl}==1));
        exptiming.respmat_sh{i_bl}=zeros(size(exptiming.respmat{i_bl}));
        exptiming.respmat_sh{i_bl}(find(exptiming.respmat{i_bl})+t.shifts)=1;
        
        t.shifts = arrayfun(@(x) subsref(strfind(exptiming.framemat{i_bl}(x-round(numel(F.flickframes)/2):end),F.flickframes(1:end/2)),...
            struct('type','()','subs',{{1}}))-round(numel(F.flickframes)/2)-1,...
            find(exptiming.blinkmat{i_bl}==1));
        exptiming.blinkmat_sh{i_bl}=zeros(size(exptiming.blinkmat{i_bl}));
        exptiming.blinkmat_sh{i_bl}(find(exptiming.blinkmat{i_bl})+t.shifts)=1;
        
        % t.tt=[exptiming.frametime{i_bl}; exptiming.framemat{i_bl};  exptiming.triggermat{i_bl}; exptiming.respmat{i_bl}; exptiming.respmat_sh{i_bl}];
    end
    
    %% write new events
    for i_bl = 1:6
        % shifted response events
        t.ev=find(exptiming.respmat_sh{i_bl});
        for i_ev = 1:numel(t.ev)
            EEG_exp = pop_editeventvals(EEG_exp,'append',{1 [] [] [] 1},...
                'changefield',{2 'epoch' i_bl},...
                'changefield',{2 'latency' exptiming.frametime{i_bl}(t.ev(i_ev))},...
                'changefield',{2 'type' str2double(F.RespTrigger{1})+5});
        end
        
        % shifted blinks
        t.ev=find(exptiming.blinkmat_sh{i_bl});
        for i_ev = 1:numel(t.ev)
            EEG_exp = pop_editeventvals(EEG_exp,'append',{1 [] [] [] 1},...
                'changefield',{2 'epoch' i_bl},...
                'changefield',{2 'latency' exptiming.frametime{i_bl}(t.ev(i_ev))},...
                'changefield',{2 'type' str2double(F.BlinkTrigger{1})+5});
        end
        
        % numel(find(cell2mat(EEG_exp.epoch(6).eventtype)==30))
        % numel(find(cell2mat(EEG_exp.epoch(6).eventtype)==35))
    end
    EEG_exp=eeg_checkset(EEG_exp);
   
    %% extract new epocs + artifact correction
    % pop_eegplot(EEG_exp,1,1,1)
            
    % chanlocs
    EEG_exp.chanlocs=pop_chanedit(EEG_exp.chanlocs,'load',{F.ChanlocFile,'filetype','besa (elp)'}); % Load Channel Locations
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    EEG_ext = pop_select(EEG_exp,'channel',65:72);
    EEG_exp = pop_select(EEG_exp,'channel',1:64);
    
    if F.CSD_flag == 1 & i_sub == 1 % calculate CSD matrix
        CSD.chanmat=ExtractMontage('C:\Users\psy05cvd\Dropbox\work\matlab\software\toolboxes\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG_exp.chanlocs.labels}');
        [CSD.G,CSD.H] = GetGH(CSD.chanmat);
    end
    
    % epoch (depending on F.TFAFlag)
    EEG_exp_resp_ep_n = pop_epoch( EEG_exp, F.RespTrigger, F.RespEpoch, 'epochinfo', 'yes');
    EEG_exp_resp_ep_sh = pop_epoch( EEG_exp, {num2str(str2double(F.RespTrigger{1})+5)}, F.RespEpoch, 'epochinfo', 'yes');
    EEG_exp_blink_ep_n = pop_epoch( EEG_exp, F.BlinkTrigger, F.EOGEpoch, 'epochinfo', 'yes');
    EEG_exp_blink_ep_sh = pop_epoch( EEG_exp, {num2str(str2double(F.BlinkTrigger{1})+5)}, F.EOGEpoch, 'epochinfo', 'yes');
    EEG_ext_resp_ep_n = pop_epoch( EEG_ext, F.RespTrigger, F.RespEpoch, 'epochinfo', 'yes');
    EEG_ext_resp_ep_sh = pop_epoch( EEG_ext, {num2str(str2double(F.RespTrigger{1})+5)}, F.RespEpoch, 'epochinfo', 'yes');
    EEG_ext_blink_ep_n = pop_epoch( EEG_ext, F.BlinkTrigger, F.EOGEpoch, 'epochinfo', 'yes');
    EEG_ext_blink_ep_sh = pop_epoch( EEG_ext, {num2str(str2double(F.BlinkTrigger{1})+5)}, F.EOGEpoch, 'epochinfo', 'yes');
    
    
    % get timing
    % get timing (adjusted)
    art.resp_time_n=[];
    for i_tr = 1:EEG_exp_resp_ep_n.trials
        t.ix=EEG_exp_resp_ep_n.epoch(i_tr).eventurevent{dsearchn(cell2mat(EEG_exp_resp_ep_n.epoch(i_tr).eventlatency)',0)};
        art.resp_time_n(i_tr)=EEG_exp_resp_ep_n.urevent(t.ix).latency./EEG_exp_resp_ep_n.srate;
    end
    art.resp_time_sh=[];
    for i_tr = 1:EEG_exp_resp_ep_sh.trials
        t.ix=EEG_exp_resp_ep_sh.epoch(i_tr).eventurevent{dsearchn(cell2mat(EEG_exp_resp_ep_sh.epoch(i_tr).eventlatency)',0)};
        art.resp_time_sh(i_tr)=EEG_exp_resp_ep_sh.urevent(t.ix).latency./EEG_exp_resp_ep_sh.srate;
    end
    art.blink_time_n=[];
    for i_tr = 1:EEG_exp_blink_ep_n.trials
        t.ix=EEG_exp_blink_ep_n.epoch(i_tr).eventurevent{dsearchn(cell2mat(EEG_exp_blink_ep_n.epoch(i_tr).eventlatency)',0)};
         art.blink_time_n(i_tr)=EEG_exp_blink_ep_n.urevent(t.ix).latency./EEG_exp_blink_ep_n.srate;
    end
    art.blink_time_sh=[];
    for i_tr = 1:EEG_exp_blink_ep_sh.trials
        t.ix=EEG_exp_blink_ep_sh.epoch(i_tr).eventurevent{dsearchn(cell2mat(EEG_exp_blink_ep_sh.epoch(i_tr).eventlatency)',0)};
        art.blink_time_sh(i_tr)=EEG_exp_blink_ep_sh.urevent(t.ix).latency./EEG_exp_blink_ep_sh.srate;
    end
    
    % see in which trials blinks are found
    art.blink_time_in_resptrials_n =[];
    art.blinknum_in_resptrials_n =[];
    for i_tr = 1:EEG_exp_resp_ep_n.trials
        tix = cell2mat(EEG_exp_resp_ep_n.epoch(i_tr).eventtype)==str2num(F.BlinkTrigger{1});
        art.blinknum_in_resptrials_n(i_tr) = sum(tix);
        art.blink_time_in_resptrials_n{i_tr,1}= cell2mat(EEG_exp_resp_ep_n.epoch(i_tr).eventlatency(tix));
    end
    
    art.blink_time_in_resptrials_sh =[];
    art.blinknum_in_resptrials_sh =[];
    for i_tr = 1:EEG_exp_resp_ep_sh.trials
        tix = cell2mat(EEG_exp_resp_ep_sh.epoch(i_tr).eventtype)==str2num(F.BlinkTrigger{1});
        art.blinknum_in_resptrials_sh(i_tr) = sum(tix);
        art.blink_time_in_resptrials_sh{i_tr,1}= cell2mat(EEG_exp_resp_ep_sh.epoch(i_tr).eventlatency(tix));
    end
    
    
    % extract blink time in each trial
    art.blink_trialtimes_n = [];
    for i_tr = 1:EEG_exp_resp_ep_n.trials
        t.idx=cellfun(@(x) x==str2double(F.BlinkTrigger),EEG_exp_resp_ep_n.epoch(i_tr).eventtype);
        art.blink_trialtimes_n = [art.blink_trialtimes_n EEG_exp_resp_ep_n.epoch(i_tr).eventlatency{t.idx}];
    end
    art.blink_trialtimes_sh = [];
    for i_tr = 1:EEG_exp_resp_ep_sh.trials
        t.idx=cellfun(@(x) x==str2double(F.BlinkTrigger),EEG_exp_resp_ep_sh.epoch(i_tr).eventtype);
        art.blink_trialtimes_sh = [art.blink_trialtimes_sh EEG_exp_resp_ep_sh.epoch(i_tr).eventlatency{t.idx}];
    end
    % figure; hist(art.blink_trialtimes_sh,30);
    
      
    % detrend
    EEG_exp_resp_ep_n = eegF_Detrend( EEG_exp_resp_ep_n);
    EEG_exp_resp_ep_sh = eegF_Detrend( EEG_exp_resp_ep_sh);
    EEG_exp_blink_ep_n = eegF_Detrend( EEG_exp_blink_ep_n);
    EEG_exp_blink_ep_sh = eegF_Detrend( EEG_exp_blink_ep_sh);
    EEG_ext_resp_ep_n = eegF_Detrend( EEG_ext_resp_ep_n);
    EEG_ext_resp_ep_sh = eegF_Detrend( EEG_ext_resp_ep_sh);
    EEG_ext_blink_ep_n = eegF_Detrend( EEG_ext_blink_ep_n);
    EEG_ext_blink_ep_sh = eegF_Detrend( EEG_ext_blink_ep_sh);
    
    %pop_eegplot(EEG_exp_resp_ep_n,1,1,1)
    %pop_eegplot(EEG_exp_resp_ep_sh,1,1,1)
    %figure; plot(EEG_ext_resp_ep_n.times,mean(EEG_ext_resp_ep_n.data,3)); legend({EEG_ext_resp_ep_n.chanlocs.labels})
    %figure; plot(EEG_ext_resp_ep_sh.times,mean(EEG_ext_resp_ep_sh.data,3)); legend({EEG_ext_resp_ep_sh.chanlocs.labels})
    %figure; plot(EEG_ext_resp_ep_n.times,std(EEG_ext_resp_ep_n.data,1,3)); legend({EEG_ext_resp_ep_n.chanlocs.labels})
    %figure; plot(EEG_ext_resp_ep_sh.times,std(EEG_ext_resp_ep_sh.data,1,3)); legend({EEG_ext_resp_ep_sh.chanlocs.labels})
    %figure; plot(EEG_ext_resp_ep_n.times,mean(EEG_ext_resp_ep_n.data(1,:,:)-EEG_ext_resp_ep_n.data(2,:,:),3))
    %figure; plot(EEG_ext_resp_ep_n.times,std(EEG_ext_resp_ep_n.data(1,:,:)-EEG_ext_resp_ep_n.data(2,:,:),1,3))
   
    %pop_eegplot(EEG_expep,1,1,1)
    %figure; pop_plottopo(EEG_exp_resp_ep_n, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    %figure; pop_plottopo(EEG_exp_resp_ep_sh, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    %figure; pop_plottopo(EEG_exp_blink_ep_n, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    %figure; pop_plottopo(EEG_exp_blink_ep_sh, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    
    %% calculate RESS component
    [RESS.exp_resp_ep_sh] = RESS_Calculate_Component(EEG_exp_resp_ep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqdist', 2, 'neighfreqwidth', 2);
    %[RESS.exp_resp_ep_sh] = RESS_Calculate_Component(EEG_exp_resp_ep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqdist', 7, 'neighfreqwidth', 5); % best for sub 15
    %[RESS.exp_resp_ep_sh] = RESS_Calculate_Component(EEG_exp_resp_ep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqdist', 7, 'neighfreqwidth', 6);
%     [RESS.exp_resp_ep_sh] = RESS_Calculate_Component(EEG_exp_resp_ep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqdist', 1, 'neighfreqwidth', 1,'diagnostics','on');


    
    
    % save componenet
    EOG.RESS.map=RESS.exp_resp_ep_sh{1}.ComponentMaps(:,RESS.exp_resp_ep_sh{1}.CompNum)./...
        max(RESS.exp_resp_ep_sh{1}.ComponentMaps(:,RESS.exp_resp_ep_sh{1}.CompNum));
    EOG.RESS.map_signchanged=RESS.exp_resp_ep_sh{1}.ComponentMaps_signchanged(:,RESS.exp_resp_ep_sh{1}.CompNum)./...
        max(RESS.exp_resp_ep_sh{1}.ComponentMaps_signchanged(:,RESS.exp_resp_ep_sh{1}.CompNum));
    EOG.RESS.SNR_ind=RESS.exp_resp_ep_sh{1}.Comp_SNR_induced;
    EOG.RESS.SNR_evo=RESS.exp_resp_ep_sh{1}.Comp_SNR_evoked;
    EOG.RESS.eigenvectors=RESS.exp_resp_ep_sh{1}.eigenvectors;
    EOG.RESS.eigenvalues=RESS.exp_resp_ep_sh{1}.eigenvalues;
    EOG.RESS.parameters=RESS.exp_resp_ep_sh{1}.parameters;
    EOG.RESS.trialnum = EEG_exp_resp_ep_sh.trials;
    
    % figure; topoplot(EOG.RESSmap,EEG_exp_resp_ep_sh.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off');
    % figure; topoplot(EOG.RESSmap_signchanged,EEG_exp_resp_ep_sh.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off');
    
    EEG_exp_resp_ep_sh_ress = RESS.exp_resp_ep_sh{1,1}.EEG_out;
        
    
    % use spatial filters of shifted data for not shifted response data
    EEG_exp_resp_ep_n_ress=EEG_exp_resp_ep_n;
     
    EEG_exp_resp_ep_n_ress = pop_select( EEG_exp_resp_ep_n_ress,'channel',1);
    EEG_exp_resp_ep_n_ress=pop_chanedit(EEG_exp_resp_ep_n_ress, 'changefield',{1 'labels' 'RESS'},...
        'changefield',{1 'theta' ''},'changefield',{1 'radius' ''},'changefield',...
        {1 'X' ''},'changefield',{1 'Y' ''},'changefield',{1 'Z' ''},'changefield',{1 'sph_theta' ''},...
        'changefield',{1 'sph_phi' ''},'changefield',{1 'sph_radius' ''});
   
    ress_ts = zeros(size(EEG_exp_resp_ep_n_ress.data));
    for i_tr=1:EEG_exp_resp_ep_n_ress.trials
        ress_ts(1,:,i_tr) = RESS.exp_resp_ep_sh{1}.eigenvectors(:,RESS.exp_resp_ep_sh{1}.CompNum)'*squeeze(EEG_exp_resp_ep_n.data(:,:,i_tr));
    end
    EEG_exp_resp_ep_n_ress.data=ress_ts;
    
    % use spatial filters for shifted blink data
    EEG_exp_blink_ep_sh_ress = EEG_exp_blink_ep_sh;
    EEG_exp_blink_ep_sh_noress = EEG_exp_blink_ep_sh;
    EEG_exp_blink_ep_n_ress = EEG_exp_blink_ep_n;
    EEG_exp_blink_ep_n_noress = EEG_exp_blink_ep_n;
     
    EEG_exp_blink_ep_sh_ress = pop_select( EEG_exp_blink_ep_sh_ress,'channel',1);
    EEG_exp_blink_ep_sh_ress = pop_chanedit(EEG_exp_blink_ep_sh_ress, 'changefield',{1 'labels' 'RESS'},...
        'changefield',{1 'theta' ''},'changefield',{1 'radius' ''},'changefield',...
        {1 'X' ''},'changefield',{1 'Y' ''},'changefield',{1 'Z' ''},'changefield',{1 'sph_theta' ''},...
        'changefield',{1 'sph_phi' ''},'changefield',{1 'sph_radius' ''});
    EEG_exp_blink_ep_n_ress = pop_select( EEG_exp_blink_ep_n_ress,'channel',1);
    EEG_exp_blink_ep_n_ress = pop_chanedit(EEG_exp_blink_ep_n_ress, 'changefield',{1 'labels' 'RESS'},...
        'changefield',{1 'theta' ''},'changefield',{1 'radius' ''},'changefield',...
        {1 'X' ''},'changefield',{1 'Y' ''},'changefield',{1 'Z' ''},'changefield',{1 'sph_theta' ''},...
        'changefield',{1 'sph_phi' ''},'changefield',{1 'sph_radius' ''});
   
    ress_ts = zeros(size(EEG_exp_blink_ep_sh_ress.data));
    for i_tr=1:EEG_exp_blink_ep_sh_ress.trials
        ress_ts(1,:,i_tr) = RESS.exp_resp_ep_sh{1}.eigenvectors(:,RESS.exp_resp_ep_sh{1}.CompNum)'*squeeze(EEG_exp_blink_ep_sh.data(:,:,i_tr));
    end
    EEG_exp_blink_ep_sh_ress.data=ress_ts;
    ress_ts = zeros(size(EEG_exp_blink_ep_n_ress.data));
    for i_tr=1:EEG_exp_blink_ep_n_ress.trials
        ress_ts(1,:,i_tr) = RESS.exp_resp_ep_sh{1}.eigenvectors(:,RESS.exp_resp_ep_sh{1}.CompNum)'*squeeze(EEG_exp_blink_ep_n.data(:,:,i_tr));
    end
    EEG_exp_blink_ep_n_ress.data=ress_ts;
    
    % noress data
    noress_ts = zeros(size(EEG_exp_blink_ep_sh_noress.data));
    for i_tr=1:EEG_exp_blink_ep_sh_noress.trials
        noress_ts(:,:,i_tr)  = RESS.exp_resp_ep_sh{1}.ComponentMaps(:,1:end-1) * RESS.exp_resp_ep_sh{1}.eigenvectors(:,1:end-1)' * EEG_exp_blink_ep_sh_noress.data(:,:,i_tr);
    end
    EEG_exp_blink_ep_sh_noress.data=noress_ts;
    
    noress_ts = zeros(size(EEG_exp_blink_ep_n_noress.data));
    for i_tr=1:EEG_exp_blink_ep_n_noress.trials
        noress_ts(:,:,i_tr)  = RESS.exp_resp_ep_sh{1}.ComponentMaps(:,1:end-1) * RESS.exp_resp_ep_sh{1}.eigenvectors(:,1:end-1)' * EEG_exp_blink_ep_n_noress.data(:,:,i_tr);
    end
    EEG_exp_blink_ep_n_noress.data=noress_ts;
    
    %% %% calculate current source density transform
    fprintf(1,'\n###\ncalculating CSD transform\n###\n')
    if F.CSD_flag == 1
        for i_tr = 1:EEG_exp_blink_ep_sh_noress.trials
            % csd of raw data
            EEG_exp_blink_ep_sh_noress.data(:,:,i_tr)= CSDTransform(EEG_exp_blink_ep_sh_noress.data(:,:,i_tr), CSD.G, CSD.H);
        end
        for i_tr = 1:EEG_exp_blink_ep_n_noress.trials
            % csd of raw data
            EEG_exp_blink_ep_n_noress.data(:,:,i_tr)= CSDTransform(EEG_exp_blink_ep_n_noress.data(:,:,i_tr), CSD.G, CSD.H);
        end
    end
    
    %% do TFA
    TFA.electrodes = EEG_exp_blink_ep_sh_noress.chanlocs;
    TFA.params.srate = EEG_exp_blink_ep_sh_noress.srate;
    
    
    %filter
    TFA.params.filter = {[] []};
    %TFA.params.filter = {0.5 []};
    %EEG_expep_n = pop_eegfiltnew(EEG_expep_n, TFA.params.filter{1}, TFA.params.filter{2}, EEG_expep_n.pnts, 0, [], 0);
    %EEG_expep_sh = pop_eegfiltnew(EEG_expep_sh, TFA.params.filter{1}, TFA.params.filter{2}, EEG_expep_sh.pnts, 0, [], 0);
    
    TFA.data_induced=nan(numel(F.TFAfreqs),EEG_exp_blink_ep_n_noress.pnts/2, numel(F.EEGChans));
    TFA.data_evoked=TFA.data_induced;
    TFA.data_RESS_induced = nan(numel(F.TFAfreqs),EEG_exp_blink_ep_n_noress.pnts/2);
    TFA.data_RESS_evoked=TFA.data_RESS_induced;
    TFA.alltrials = EEG_exp_blink_ep_n_noress.trials;
    
    
    fprintf(1,'calculating %1.0f TFAs:\n',numel(F.TFAfreqs))
    TFA.type = 'gabor';
    TFA.f = F.TFAfreqs;
    %     TFA.params.gabor_FWHM_freq = 0.5;
    TFA.params.gabor_FWHM_freq = 1;
    for i_freq = 1:numel(TFA.f)
        fprintf(1,'%1.0f of %1.0f\n',i_freq, numel(F.TFAfreqs))
        % induced without RESS component
        EEG_Gabor_gab = eegF_Gabor(EEG_exp_blink_ep_sh_noress, TFA.f(i_freq),TFA.params.gabor_FWHM_freq);
        EEG_Gabor_gab=pop_select( EEG_Gabor_gab,'time',F.EOGEpoch);
        TFA.data_induced(i_freq,:,:)=imresize(mean(EEG_Gabor_gab.data,3),[EEG_Gabor_gab.nbchan,EEG_Gabor_gab.pnts/2])';
        t.times = EEG_Gabor_gab.times([1 end]);
        
        % evoked without RESS component
        EEG_t = EEG_exp_blink_ep_sh_noress; EEG_t=pop_select(EEG_t,'trial',1); EEG_t.data=mean(EEG_exp_blink_ep_sh_noress.data,3);
        EEG_Gabor_gab = eegF_Gabor(EEG_t, TFA.f(i_freq),TFA.params.gabor_FWHM_freq);
        TFA.data_evoked(i_freq,:,:)=imresize(...
            EEG_Gabor_gab.data(:,dsearchn(EEG_Gabor_gab.times',t.times(1)):dsearchn(EEG_Gabor_gab.times',t.times(2))),...
            [EEG_Gabor_gab.nbchan,numel(dsearchn(EEG_Gabor_gab.times',t.times(1)):dsearchn(EEG_Gabor_gab.times',t.times(2)))/2])';
        
        % evoked RESS component
        EEG_t = EEG_exp_blink_ep_sh_ress; EEG_t=pop_select(EEG_t,'trial',1); EEG_t.data=mean(EEG_exp_blink_ep_sh_ress.data,3);
        EEG_Gabor_gab = eegF_Gabor(EEG_t, TFA.f(i_freq),TFA.params.gabor_FWHM_freq);
        TFA.data_RESS_evoked(i_freq,:,:)=imresize(...
            EEG_Gabor_gab.data(:,dsearchn(EEG_Gabor_gab.times',t.times(1)):dsearchn(EEG_Gabor_gab.times',t.times(2))),...
            [EEG_Gabor_gab.nbchan,numel(dsearchn(EEG_Gabor_gab.times',t.times(1)):dsearchn(EEG_Gabor_gab.times',t.times(2)))/2])';
        
        % induced RESS component
        EEG_Gabor_gab = eegF_Gabor(EEG_exp_blink_ep_sh_ress, TFA.f(i_freq),TFA.params.gabor_FWHM_freq);
        EEG_Gabor_gab=pop_select( EEG_Gabor_gab,'time',F.EOGEpoch);
        TFA.data_RESS_induced(i_freq,:,:)=imresize(mean(EEG_Gabor_gab.data,3),[EEG_Gabor_gab.nbchan,EEG_Gabor_gab.pnts/2])';
    end
        
            
            % plotting for checking
%             pl.elec2plot = {'Oz'};
% %             pl.elec2plot = {'C3';'CP3'};
%             pl.base=[-3500 3000]; pl.base_i=dsearchn(TFA.t',pl.base');
%             pl.t2pl = [-4000 4000]; pl.t2pl_i = dsearchn(TFA.t',pl.t2pl');
%             pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
%             
%             pl.data = TFA.data_RESS_induced(:,pl.t2pl_i(1):pl.t2pl_i(2),:); pl.ti1 = 'RESS'; pl.ti22 = 'raw';  pl.ti3 = 'induced'; pl.clims = [0 max(pl.data(:))];
%             pl.data = TFA.data_RESS_evoked(:,pl.t2pl_i(1):pl.t2pl_i(2),:); pl.ti1 = 'RESS'; pl.ti2 = 'raw'; pl.ti3 = 'evoked'; pl.clims = [0 max(pl.data(:))];
%             pl.data = mean(TFA.data_induced(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i),3); pl.ti1 = pl.elec2plot; pl.ti2 = 'raw';  pl.ti3 = 'induced'; pl.clims = [0 max(pl.data(:))];
%             pl.data = mean(TFA.data_evoked(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i),3); pl.ti1 = pl.elec2plot; pl.ti2 = 'raw'; pl.ti3 = 'evoked'; pl.clims = [0 max(pl.data(:))];
%             
%             pl.data = bsxfun(@minus, TFA.data_RESS_induced(:,pl.t2pl_i(1):pl.t2pl_i(2),:), mean(TFA.data_RESS_induced(:,pl.base_i(1):pl.base_i(2)),2));
%             pl.ti1 = 'RESS'; pl.ti22 = 'baseline corrected';  pl.ti3 = 'induced'; pl.clims = [-1 1]*max(abs(pl.data(:)));
%             pl.data = bsxfun(@minus, TFA.data_RESS_evoked(:,pl.t2pl_i(1):pl.t2pl_i(2),:), mean(TFA.data_RESS_evoked(:,pl.base_i(1):pl.base_i(2)),2));
%             pl.ti1 = 'RESS'; pl.ti22 = 'baseline corrected';  pl.ti3 = 'evoked'; pl.clims = [-1 1]*max(abs(pl.data(:)));
%             pl.data = mean(bsxfun(@minus, TFA.data_induced(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i), mean(TFA.data_induced(:,pl.base_i(1):pl.base_i(2),pl.elec2plot_i),2)),3);
%             pl.ti1 = pl.elec2plot; pl.ti22 = 'baseline corrected';  pl.ti3 = 'induced'; pl.clims = [-1 1]*max(abs(pl.data(:)));
%             pl.data = mean(bsxfun(@minus, TFA.data_evoked(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i), mean(TFA.data_evoked(:,pl.base_i(1):pl.base_i(2),pl.elec2plot_i),2)),3);
%             pl.ti1 = pl.elec2plot; pl.ti22 = 'baseline corrected';  pl.ti3 = 'evoked'; pl.clims = [-1 1]*max(abs(pl.data(:)));
%       
%             figure;
%             imagesc(TFA.t(pl.t2pl_i(1):pl.t2pl_i(2)),TFA.f,pl.data,pl.clims)
%             set(gca,'YDir','normal')
%             title(sprintf('%s %s tfa %s', pl.ti2, pl.ti3, vararg2str(pl.ti1)), 'FontSize',8)
%             colorbar
%             xlabel('time in ms')
%             ylabel('frequency in Hz')
%             hline(14.16667,'m')
%             set(gca,'FontSize',8)
    
    fprintf(1,'...done\n')
    TFA.t=imresize(EEG_Gabor_gab.times,[1,EEG_Gabor_gab.pnts/2]);
    
    %% additional analysis steps
    pl.xdata = EEG_exp_blink_ep_sh_ress.times;
    pl.ydata1 = mean(EEG_ext_blink_ep_sh.data(1,:,:)-EEG_ext_blink_ep_sh.data(2,:,:),3);
    pl.ydata2 = mean(EEG_exp_blink_ep_sh_ress.data,3);
    [t.tfadata,t.t,t.f,t.D] = mwavelTFA(squeeze(EEG_exp_blink_ep_sh_ress.data),85/6,EEG_exp_blink_ep_sh_ress.srate,(85/6)/3,1,3.5);
    pl.ydata3 = mean(t.tfadata,3);
    [t.tfadata,t.t,t.f,t.D] = mwavelTFA(squeeze(mean(EEG_exp_blink_ep_sh_ress.data,3)),85/6,EEG_exp_blink_ep_sh_ress.srate,(85/6)/3,1,3.5);
    pl.ydata4 = t.tfadata;
    fig1 = figure; set(fig1, 'WindowStyle', 'docked')
    subplot(2,1,1)
    [h.AX,h.line1,h.line2] = plotyy(pl.xdata,pl.ydata1,repmat(pl.xdata',1,3),[pl.ydata2' pl.ydata3' pl.ydata4']);
    grid on
    set(h.AX(1),'ylim',[-1 1]*max(abs(get(h.AX(1),'ylim'))),'YTickMode','auto')
    t.max = max(max(max(abs([pl.ydata2' pl.ydata3' pl.ydata4']))));
    %t.sc = ceil(t.max*(10^(fix(1/t.max)-1)))/10^(fix(1/t.max)-1);
    if t.max<1
        t.sc = ceil(t.max*10^numel(num2str(fix(1/t.max))))/(10^numel(num2str(fix(1/t.max))));
    else
        t.sc = ceil(t.max*1)/1;
    end
    
    set(h.AX(2),'ylim',[-1 1]*t.sc,'YTickMode','auto')
    title('RESS')
%     legend({'VEOG';'SSVEP signal';'induced amp';'evoked amp'},'Orientation','horizontal','Location','SouthOutside')
    
    pl.xdata = EEG_exp_blink_ep_n_ress.times;
    pl.ydata1 = mean(EEG_ext_blink_ep_n.data(1,:,:)-EEG_ext_blink_ep_n.data(2,:,:),3);
    pl.ydata2 = mean(EEG_exp_blink_ep_n_ress.data,3);
    [t.tfadata,t.t,t.f,t.D] = mwavelTFA(squeeze(EEG_exp_blink_ep_n_ress.data),85/6,EEG_exp_blink_ep_n_ress.srate,(85/6)/3,1,3.5);
    pl.ydata3 = mean(t.tfadata,3);
    [t.tfadata,t.t,t.f,t.D] = mwavelTFA(squeeze(mean(EEG_exp_blink_ep_n_ress.data,3)),85/6,EEG_exp_blink_ep_n_ress.srate,(85/6)/3,1,3.5);
    pl.ydata4 = t.tfadata;
    subplot(2,1,2)
    [h.AX,h.line1,h.line2] = plotyy(pl.xdata,pl.ydata1,repmat(pl.xdata',1,3),[pl.ydata2' pl.ydata3' pl.ydata4']);
    grid on
    set(h.AX(1),'ylim',[-1 1]*max(abs(get(h.AX(1),'ylim'))),'YTickMode','auto')
    t.max = max(max(max(abs([pl.ydata2' pl.ydata3' pl.ydata4']))));
    if t.max<1
        t.sc = ceil(t.max*10^numel(num2str(fix(1/t.max))))/(10^numel(num2str(fix(1/t.max))));
    else
        t.sc = ceil(t.max*1)/1;
    end
    set(h.AX(2),'ylim',[-1 1]*t.sc,'YTickMode','auto')
    
    legend({'VEOG';'SSVEP signal';'induced amp';'evoked amp'},'Orientation','horizontal','Location','SouthOutside')
    title('NoRess')
        
    %% save EOG files
    EOG.EEG_ext_blink_ep_n = EEG_ext_blink_ep_n;
    EOG.EEG_ext_blink_ep_sh = EEG_ext_blink_ep_sh;
    EOG.EEG_ext_resp_ep_n = EEG_ext_resp_ep_n;
    EOG.EEG_ext_resp_ep_sh = EEG_ext_resp_ep_sh;
    EOG.EEG_exp_blink_ep_n_ress = EEG_exp_blink_ep_n_ress;
    EOG.EEG_exp_blink_ep_sh_ress = EEG_exp_blink_ep_sh_ress;
    EOG.TFA = TFA;
    
    %% save
    EOG.art=art;
    EOG.savetime=datestr(now);
    fprintf(1,'|| saving file ||  %s\\VP%02.0f_exp_EOG.mat ||\n', F.PathOut,F.Subjects2Use(i_sub)')
    save(sprintf('%s\\VP%02.0f_exp_EOG.mat', F.PathOut,F.Subjects2Use(i_sub)), 'EOG')   
    
end
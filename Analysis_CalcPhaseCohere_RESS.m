%% script in order to calculate time frequency analysis of main experiment



clearvars

%% parameters
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\RAW';
F.PathOut               = 'D:\work\data\SSVEP_volmov\EEG\PHASECOHERE_RESS_noblinks_CSD';
F.PathInEOG             = 'D:\work\data\SSVEP_volmov\EEG\EOG_Data\';
% F.PathInRESS            = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD2';
F.PathInRESS            = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD4'; % gabor of 1Hz FWHM [no blinks -3.5 3.5]
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27];
F.Subjects2Use          = [3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27];
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.VEOGChans             = [];
F.HEOGChans             = [];
F.ChanlocFile           = 'C:\Users\psy05cvd\Dropbox\work\matlab\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.RespTrigger           = {'30'};
F.RespEpoch             = [-5.5 5.5]; % in s
F.SSVEPTrigger          = {'15'};

F.ExpEpochTrigger       = {'11','12'};
F.eptime                = 360; % time of each block in s
F.Blocks                = 6; % number of blocks

F.FrameRate             = 85;
F.TriggerRate           = 6; % SSVEP trigger frequency (e.g. 6 = every 6s)
F.flickframes           = [1 1 1 0 0 0]; % on-off frames for SSVEP

%F.TFAfreqs              = [5:(1/6):40];
F.TFAfreqs              = [4-(1/3):0.25:45];
F.TFAFlag               = [2]; % 1=wavelet; 2=gabor
F.TFA_FWHM              = [1];
F.CSD_flag              = 1; % 0 = no; 1 = yes

%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp.bdf ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    if F.Subjects2Use(i_sub) ~= 24
        EEG = pop_readbdf(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)), [] ,F.EMGChans(2)+1,[]);
        % pop_eegplot(EEG,1,1,1);
    else
        
        % troubleshooting VP24
        EEG = pop_loadset('filename','VP24_exp.set','filepath',F.PathIn);
    end
    
    % pop_eegplot(EEG,1,1,1);
    %numel(find(cell2mat({EEG.event.type})==15))
    
    % load EOG data
    EOG = open(sprintf('%s\\VP%02.0f_exp_EOG.mat',F.PathInEOG,F.Subjects2Use(i_sub)));
    
    % rereference
    EEG = pop_reref( EEG, [],'exclude',F.EEGChans(end)+1:EEG.nbchan ); % average
    % pop_eegplot(EEG,1,1,1,1)
    
%     EEG = pop_resample(EEG, 256);
    
    %% TFA power around button press
    % induced and evoked...requires shift of button press to SSVEP cycle beginning
    
    % index block length
    t.i1=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{1});
    t.i2=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{2});
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
    
    %% loop for blocks
    % create time vectors for visual stimulation
    t.bl=1;
    for i_bl = 1:F.Blocks
        % trigger times for
        t.t_trig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.SSVEPTrigger{1})}}); % SSVEP trigger
        t.t_btrig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.RespTrigger{1})}}); % behavior trigger
        
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
    
    
    % discard trials with blinks around movement period [-3 +3] s
    t.rejindex_n = false(1,EEG_expep_n.trials);
    t.rejindex_sh = false(1,EEG_expep_sh.trials);
    % index trials in certain time interval
    art.win_blinks = [-3500 3500];
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
    
    
    
    % run SCADS to statistically detetct artifacts for unshifted trials
%     [EEG_expep_n art.SCADS_Trials2Del_n art.SCADS_SumData_n] = SCADS_Pass1(EEG_expep_n,F.RespEpoch,[6 6 15],1);
    % discard trials indexed as containing artifacts
%     t.rejindex = EEG_expep_n.reject.rejmanual;
%     EEG_expep_n = pop_rejepoch(EEG_expep_n, EEG_expep_n.reject.rejmanual, 0);
    
    % run SCADS to statistically detetct artifacts for shifted trials
%     [EEG_expep_sh art.SCADS_Trials2Del_sh art.SCADS_SumData_sh] = SCADS_Pass1(EEG_expep_sh,F.RespEpoch,[6 6 15],1);
    % discard trials indexed as containing artifacts
%     EEG_expep_sh = pop_rejepoch(EEG_expep_sh, EEG_expep_sh.reject.rejmanual, 0);
    
    % reject trials for external signals
%     EEG_ext_rej = pop_rejepoch(EEG_ext, t.rejindex, 0);
    EEG_ext_rej = EEG_ext;
    %figure; plot(EEG_ext_rej.times,mean(EEG_ext_rej.data,3)); legend({EEG_ext_rej.chanlocs.labels})
    %figure; plot(EEG_ext_rej.times,std(EEG_ext_rej.data,1,3)); legend({EEG_ext_rej.chanlocs.labels})
   
    % detrend
    EEG_expep_n = eegF_Detrend(EEG_expep_n);
    EEG_expep_sh = eegF_Detrend(EEG_expep_sh);
    EEG_ext = eegF_Detrend(EEG_ext);
    %pop_eegplot(EEG_expep_n,1,1,1)
    %pop_eegplot(EEG_expep_sh,1,1,1)
    %figure; plot(EEG_ext.times,mean(EEG_ext.data,3)); legend({EEG_ext.chanlocs.labels})
    %figure; plot(EEG_ext.times,std(EEG_ext.data,1,3)); legend({EEG_ext.chanlocs.labels})
   
    %pop_eegplot(EEG_expep,1,1,1)
    %figure; pop_plottopo(EEG_expep, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    %figure; pop_plottopo(EEG_expep_sh, F.EEGChans , 'BDF file epochs', 0, 'ydir',-1);
    
    %% calculate RESS component
    %     [RESS.expep_sh] = RESS_Calculate_Component(EEG_expep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqwidth', 0.5,'diagnostics', 'on');
    %     [RESS.expep_sh] = RESS_Calculate_Component(EEG_expep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqdist', 7, 'neighfreqwidth', 7, 'diagnostics','on');
    %     [RESS.expep_sh] = RESS_Calculate_Component(EEG_expep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqdist', 2, 'neighfreqwidth', 2, 'diagnostics','on');
    %     [RESS.expep_sh] = RESS_Calculate_Component(EEG_expep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqdist', 1, 'neighfreqwidth', 1, 'diagnostics','on');
    
    [RESS.expep_sh] = RESS_Calculate_Component(EEG_expep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqdist', 2, 'neighfreqwidth', 2, 'diagnostics','on');
    
    % save output
    if ~isdir(sprintf('%s/RESS_diagnostics',F.PathOut))
        mkdir(sprintf('%s/RESS_diagnostics',F.PathOut))
    end
    figHandles = findobj('Type', 'figure');
    for i_fig = 1:2
        set(0, 'currentfigure', figHandles(end-2+i_fig));
        set(figHandles(end-2+i_fig),'Position',[100 100 1200 675],'PaperPositionMode','auto')
        print(fullfile(F.PathOut,'RESS_diagnostics',sprintf('RESS_diagnostics_sub_%02.0f_%1.0f',F.Subjects2Use(i_sub),i_fig)),...
            '-djpeg','-r300')
        saveas(figHandles(end-2+i_fig),fullfile(F.PathOut,'RESS_diagnostics',sprintf('RESS_diagnostics_sub_%02.0f_%1.0f',F.Subjects2Use(i_sub),i_fig)),'fig')
    end
    close(figHandles(end-2+[1 2]))

    
    % save componenet
    Cohere.RESS.map=RESS.expep_sh{1}.ComponentMaps(:,RESS.expep_sh{1}.CompNum)./...
        max(RESS.expep_sh{1}.ComponentMaps(:,RESS.expep_sh{1}.CompNum));
    Cohere.RESS.map_signchanged=RESS.expep_sh{1}.ComponentMaps_signchanged(:,RESS.expep_sh{1}.CompNum)./...
        max(RESS.expep_sh{1}.ComponentMaps_signchanged(:,RESS.expep_sh{1}.CompNum));
    Cohere.RESS.ComponentMaps = RESS.expep_sh{1}.ComponentMaps;
    Cohere.RESS.SNR_ind=RESS.expep_sh{1}.Comp_SNR_induced;
    Cohere.RESS.SNR_evo=RESS.expep_sh{1}.Comp_SNR_evoked;
    Cohere.RESS.Comp_selected=RESS.expep_sh{1}.CompNum;
    Cohere.RESS.eigenvectors=RESS.expep_sh{1}.eigenvectors;
    Cohere.RESS.eigenvalues=RESS.expep_sh{1}.eigenvalues;
    Cohere.RESS.parameters=RESS.expep_sh{1}.parameters;
    Cohere.RESS.trialnum = EEG_expep_sh.trials;
    
    
    % figure; topoplot(TFA.RESSmap,EEG_expep_sh.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off');
    % figure; topoplot(TFA.RESSmap_signchanged,EEG_expep_sh.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off');
    
    EEG_expep_sh_ress = RESS.expep_sh{1,1}.EEG_out;
    EEG_expep_sh_noress = RESS.expep_sh{1,1}.EEG_out_noRESS;
    
            % plot fft results
%     figure; subplot(2,1,1); plot(EEG_expep_sh_ress.times, mean(EEG_expep_sh_ress.data(:,:,:),3)); title('RESS')
%     subplot(2,1,2); plot(EEG_expep_sh_ress.times, mean(EEG_expep_sh.data(29,:,:),3)); title('Oz')
    t.fft_1=abs(fft(squeeze(mean(EEG_expep_sh_ress.data,3)),2560,2))*2/size(EEG_expep_sh_ress.data,2);
    t.fft_2=abs(fft(squeeze(mean(EEG_expep_sh.data(:,:,:),3)),2560,2))*2/size(EEG_expep_sh.data,2);
    t.fft_freqs = ((0:size(t.fft_1,2)-1)/size(t.fft_1,2)) * EEG_expep_sh_ress.srate;
    [t.max, t.sort]=max(t.fft_2(:,dsearchn(t.fft_freqs', 85/6))); % index best electrode
    t.chan = findchannellabel(t.sort , 0);
%     figure; subplot(2,1,1); plot(t.fft_freqs,t.fft_1); xlim([0 60]);title('RESS')
%     subplot(2,1,2); plot(t.fft_freqs,t.fft_2(t.sort,:)); xlim([0 60]);title(sprintf(t.chan.StandardLabel{1}))
    
    % use spatial filters of shifted data for not shifted data
    EEG_expep_n_ress=EEG_expep_n;
    EEG_expep_n_noress=EEG_expep_n;
    
    EEG_expep_n_ress = pop_select( EEG_expep_n_ress,'channel',1);
    EEG_expep_n_ress=pop_chanedit(EEG_expep_n_ress, 'changefield',{1 'labels' 'RESS'},...
        'changefield',{1 'theta' ''},'changefield',{1 'radius' ''},'changefield',...
        {1 'X' ''},'changefield',{1 'Y' ''},'changefield',{1 'Z' ''},'changefield',{1 'sph_theta' ''},...
        'changefield',{1 'sph_phi' ''},'changefield',{1 'sph_radius' ''});
    
    % ress data
    ress_ts = zeros(size(EEG_expep_n_ress.data));
    for i_tr=1:EEG_expep_n_ress.trials
        ress_ts(1,:,i_tr) = RESS.expep_sh{1}.eigenvectors(:,RESS.expep_sh{1}.CompNum)'*squeeze(EEG_expep_n.data(:,:,i_tr));
    end
    EEG_expep_n_ress.data=ress_ts;
    
    % noress data
    noress_ts = zeros(size(EEG_expep_n_noress.data));
    for i_tr=1:EEG_expep_n_noress.trials
        noress_ts(:,:,i_tr)  = RESS.expep_sh{1}.ComponentMaps(:,1:end-1) * RESS.expep_sh{1}.eigenvectors(:,1:end-1)' * EEG_expep_n_noress.data(:,:,i_tr);
    end
    EEG_expep_n_noress.data=noress_ts;
    
    EEG_surrogate =  pop_select( EEG_expep_n_noress,'time',F.RespEpoch);
    
    %% reject trials with artifacts
    % reject trials
    EEG_expep_n = pop_rejepoch(EEG_expep_n, t.rejindex_n, 0);
    EEG_expep_n_noress = pop_rejepoch(EEG_expep_n_noress, t.rejindex_n, 0);
    EEG_expep_n_ress = pop_rejepoch(EEG_expep_n_ress, t.rejindex_n, 0);
    EEG_expep_sh = pop_rejepoch(EEG_expep_sh, t.rejindex_sh, 0);
    EEG_expep_sh_noress = pop_rejepoch(EEG_expep_sh_noress, t.rejindex_sh, 0);
    EEG_expep_sh_ress = pop_rejepoch(EEG_expep_sh_ress, t.rejindex_sh, 0);
    EEG_ext = pop_rejepoch(EEG_ext, t.rejindex_n, 0);
    
    %% calculate current source density transform
    fprintf(1,'\n###\ncalculating CSD transform\n###\n')
    if F.CSD_flag == 1
        for i_tr = 1:EEG_expep_n.trials
            % csd of raw data
            EEG_expep_n.data(:,:,i_tr)= CSDTransform(EEG_expep_n.data(:,:,i_tr), CSD.G, CSD.H);
        end
        for i_tr = 1:EEG_expep_n_noress.trials
            % csd of raw data
            EEG_expep_n_noress.data(:,:,i_tr)= CSDTransform(EEG_expep_n_noress.data(:,:,i_tr), CSD.G, CSD.H);
        end
        for i_tr = 1:EEG_expep_sh.trials
            % csd of raw data
            EEG_expep_sh.data(:,:,i_tr)= CSDTransform(EEG_expep_sh.data(:,:,i_tr), CSD.G, CSD.H);
        end
        for i_tr = 1:EEG_expep_sh_noress.trials
            % csd of raw data
            EEG_expep_sh_noress.data(:,:,i_tr)= CSDTransform(EEG_expep_sh_noress.data(:,:,i_tr), CSD.G, CSD.H);
        end
    end
    
    
    %% calculate phase coherence
    % EEG2 = EEG; EEG = EEG_expep_n_noress; EEG = EEG_expep_sh_ress;
  
%     figure; pop_newtimef( EEG_expep_sh_ress, 1, 1, [-6000  5996], [5         0.5] , ...
%         'caption', 'RESS sh ress', 'baseline',[-5000 -4000], 'freqs', [0 40], 'plotphase', 'off', 'ntimesout', 400, 'padratio', 4);
%     figure; pop_newtimef( EEG_expep_n_ress, 1, 1, [-6000  5996], [5         0.5] , ...
%         'caption', 'RESS nosh ress', 'baseline',[-5000 -4000], 'freqs', [0 40], 'plotphase', 'off', 'ntimesout', 400, 'padratio', 4);
%     figure; pop_newtimef( EEG_expep_sh_noress, 1, 13, [-6000  5996], [5         0.5] , ...
%         'topovec', 13, 'elocs', EEG_expep_sh_noress.chanlocs, 'chaninfo', EEG_expep_sh_noress.chaninfo, 'caption', 'C3 sh noress', ...
%         'baseline',[-5000 -4000], 'freqs', [0 40], 'plotphase', 'off', 'ntimesout', 400, 'padratio', 4);
%     figure; pop_newtimef( EEG_expep_sh_noress, 1, 29, [-6000  5996], [5         0.5] , ...
%         'topovec', 29, 'elocs', EEG_expep_sh_noress.chanlocs, 'chaninfo', EEG_expep_sh_noress.chaninfo, 'caption', 'Oz sh noress', ...
%         'baseline',[-5000 -4000], 'freqs', [0 40], 'plotphase', 'off', 'ntimesout', 400, 'padratio', 4);
%     figure; pop_newtimef( EEG_expep_n_noress, 1, 13, [-6000  5996], [5         0.5] , ...
%         'topovec', 13, 'elocs', EEG_expep_n_noress.chanlocs, 'chaninfo', EEG_expep_n_noress.chaninfo, 'caption', 'C3 nosh noress', ...
%         'baseline',[-5000 -4000], 'freqs', [0 40], 'plotphase', 'off', 'ntimesout', 400, 'padratio', 4);
%     figure; pop_newtimef( EEG_expep_n_noress, 1, 29, [-6000  5996], [5         0.5] , ...
%         'topovec', 29, 'elocs', EEG_expep_n_noress.chanlocs, 'chaninfo', EEG_expep_n_noress.chaninfo, 'caption', 'Oz nosh noress', ...
%         'baseline',[-5000 -4000], 'freqs', [0 40], 'plotphase', 'off', 'ntimesout', 400, 'padratio', 4);
%     
%     figure; pop_newtimef( EEG_expep_n_noress, 1, 13, [-6000  5996], [0] , 'topovec', 13, 'elocs', ...
%         EEG_expep_n_noress.chanlocs, 'chaninfo', EEG_expep_n_noress.chaninfo, 'caption', 'C3', 'baseline',[-5000 -4500], 'freqs', [0 40], ...
%         'plotphase', 'off', 'ntimesout', 400, 'padratio', 4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %eeglab approach for each channel
    % FFT approach input:
%     t.data = EEG_expep_n_noress.data(13,:,:);
%     t.data = EEG_expep_sh_ress.data(1,:,:);
%     fig1=figure; [t.ersp, t.itc, t.powbase, t.times, t.freqs, t.erspboot, t.itcboot, t.itcphase]=...
%         newtimef(t.data(:)', EEG_expep_n_noress.pnts, [EEG_expep_n_noress.times(1) EEG_expep_n_noress.times(end)], EEG_expep_n_noress.srate, [0] ,...
%         'topovec',13,'elocs',EEG_expep_n_noress.chanlocs,'caption','C3','baseline',[-4500 -4000] ,'freqs',[0 40] ,...
%         'plotphase','off','ntimesout',400,'padratio',8,'winsize', EEG_expep_sh_ress.srate);

    %loop for electrodes
%     Cohere.eeglab.params.freqs = [0 40];
%     Cohere.eeglab.params.baseline = [-4500 -4000];
%     Cohere.eeglab.params.padratio = 8;
%     

    
    
    
    % approach using eeglab toolbox and included functions
    % either standard window (EEG_expep_sh_ress.srate*2) or 'winsize', EEG_expep_sh_ress.srate,
%     fig1=figure; [Cohere.ersp, Cohere.itc, Cohere.powbase, Cohere.times, Cohere.freqs, Cohere.erspboot, Cohere.itcboot, Cohere.itcphase]  ...
%         = timef(EEG_expep_sh_ress.data(:)', EEG_expep_sh_ress.pnts, [EEG_expep_sh_ress.xmin EEG_expep_sh_ress.xmax]*1000,...
%         EEG_expep_sh_ress.srate, 20, 'timesout', EEG_expep_sh_ress.pnts,'winsize', EEG_expep_sh_ress.srate*4, 'type', {'phasecoher'},'padratio',32, 'maxfreq', 18,'baseline',-5000);
%     close(fig1)
%     %pl.data = squeeze(EEG_expep_sh.data(29,:,:));
%     %figure; [Cohere.ersp, Cohere.itc, Cohere.powbase, Cohere.times, Cohere.freqs, Cohere.erspboot, Cohere.itcboot, Cohere.itcphase]  ...
%     %    = timef(pl.data(:)', EEG_expep_sh_ress.pnts, [EEG_expep_sh_ress.xmin EEG_expep_sh_ress.xmax]*1000,...
%     %    EEG_expep_sh_ress.srate, 0, 'timesout', EEG_expep_sh_ress.pnts,'winsize', EEG_expep_sh_ress.srate, 'type', {'phasecoher'},'padratio',8, 'maxfreq', 18,'baseline',-5000);
%                            
%     tidx = dsearchn(Cohere.freqs',85/6);
% %     figure; plot(Cohere.times,Cohere.ersp(tidx,:)); title('SSVEP ersp')
% %     figure; plot(Cohere.times,Cohere.itc(tidx,:)); title('SSVEP ITC')
% %     figure; plot(Cohere.times,Cohere.itcphase(tidx,:)); title('SSVEP itc-phase')
%     
%     % this is working now!
%     Cohere.SSVEP_itcphasediff=[];
%     for i_el = 1:numel(Cohere.itc(tidx,:))-1
%         Cohere.SSVEP_itcphasediff(i_el)=rad2deg(angdiff2(Cohere.itcphase(tidx,i_el),Cohere.itcphase(tidx,i_el+1)));
%     end
%     figure; rose(diff(Cohere.itcphase(tidx,:)),10000); title('SSVEP itc-phase difference');
%     figure; plot(Cohere.times(2:end-1),Cohere.SSVEP_itcphasediff(1:end-1)); title('SSVEP abs itc-phase')
%     figure; imagesc(Cohere.times,Cohere.freqs,Cohere.ersp);colormap('jet');set(gca,'ydir','normal')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % gabor intermedley
    EEG_expep_sh_chan = pop_select(EEG_expep_sh,'channel',t.chan.StandardLabel);
    EEG_expep_sh_ress_m = EEG_expep_sh_ress; EEG_expep_sh_ress_m.data(:,:,1)=mean(EEG_expep_sh_ress.data,3);
    EEG_expep_sh_chan_m = EEG_expep_sh_chan; EEG_expep_sh_chan_m.data(:,:,1)=mean(EEG_expep_sh_chan.data,3);
    [EEG_expep_sh_ress_gabor]=eegF_Gabor(EEG_expep_sh_ress,85/6,F.TFA_FWHM,'complex');
    [EEG_expep_sh_chan_gabor]=eegF_Gabor(EEG_expep_sh_chan,85/6,F.TFA_FWHM,'complex');
    [EEG_expep_sh_ress_m_gabor]=eegF_Gabor(EEG_expep_sh_ress_m,85/6,F.TFA_FWHM,'complex');
    [EEG_expep_sh_chan_m_gabor]=eegF_Gabor(EEG_expep_sh_chan_m,85/6,F.TFA_FWHM,'complex');
    t.data = squeeze(EEG_expep_sh_ress_gabor.data);
    t.data(:,:,2) = squeeze(EEG_expep_sh_chan_gabor.data);
    
    ttf.itpc_comp=squeeze(exp(1i*angle(t.data)));
        
    ttf.itpc         = squeeze(abs(mean(ttf.itpc_comp,2)));
    ttf.prefAngle    = angle(ttf.itpc_comp);
    ttf.ersp_ind     = squeeze(mean(abs(t.data),2));
    ttf.ersp_evo     = [abs(EEG_expep_sh_ress_m_gabor.data(:,:,1))',abs(EEG_expep_sh_chan_m_gabor.data(:,:,1))'];
    
    ttf.prefAngle_diff = nan(numel(EEG_expep_sh_ress.times)-1,EEG_expep_sh_ress_gabor.trials,2);
    for i_tmp=1:numel(EEG_expep_sh_ress.times)-1
        ttf.prefAngle_diff(i_tmp,:,:)=rad2deg(angdiff2(angle(ttf.itpc_comp(i_tmp,:,:)), angle(ttf.itpc_comp(i_tmp+1,:,:))));
    end
    
    % plotting for checking
%     figure; 
%     plot(EEG_expep_sh_ress_gabor.times,ttf.itpc); legend({'RESS';'best electrode'})
%     figure; 
%     plot(EEG_expep_sh_ress_gabor.times,ttf.prefAngle(:,:,1));
%     figure; 
%     plot(EEG_expep_sh_ress_gabor.times,ttf.ersp_ind); hold on
%     plot(EEG_expep_sh_ress_gabor.times,ttf.ersp_evo);
%     legend({'RESS ind';'best electrode ind';'RESS evo';'best electrode evo'})
%     figure;
%     h = histogram(ttf.prefAngle(1,:,1),[-pi:(2*pi/100):pi]); xlim([-pi pi])
%     
%     t.Ndata = nan(numel([-pi:(2*pi/100):pi])-1,size(ttf.prefAngle,1));
%     t.tdata=angle(exp(1i*(angle(ttf.itpc_comp(:,:,1))-repmat(angle(mean(ttf.itpc_comp(:,:,1),2)),[1,size(ttf.itpc_comp,2),1]))));
%     for i_t = 1:size(ttf.prefAngle,1)
%         [t.Ndata(:,i_t),t.edges] = histcounts(t.tdata(i_t,:),[-pi:(2*pi/100):pi]);
%     end
%     t.plx=EEG_expep_sh_ress_gabor.times;
%     t.ply=t.edges(1:end-1)+(diff(t.edges)./2);
%     figure; imagesc(t.plx,rad2deg(t.ply),t.Ndata,[0 max(max(t.Ndata))]);set(gca,'ydir','normal'); ylim([-50 50])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % contrast hilbert transform and gabor and wavelet from mike cohen
    % gabor
%     [EEG_expep_sh_ress_gabor]=eegF_Gabor(EEG_expep_sh_ress,85/6,0.5,'complex');
%     [EEG_expep_sh_chan_gabor]=eegF_Gabor(EEG_expep_sh_chan,85/6,0.5,'complex');
%     t.data = squeeze(EEG_expep_sh_ress_gabor.data);
%     t.data(:,:,2) = squeeze(EEG_expep_sh_chan_gabor.data);
%     % hilbert
%     EEG_expep_sh_ress_f = pop_eegfiltnew(EEG_expep_sh_ress,(85/6)-0.5,(85/6)+0.5, 32*EEG_expep_sh_ress.srate, 0, [], 0);
%     EEG_expep_sh_chan = pop_select(EEG_expep_sh,'channel',t.chan.StandardLabel);
%     EEG_expep_sh_chan_f = pop_eegfiltnew(EEG_expep_sh_chan,(85/6)-0.5,(85/6)+0.5, 32*EEG_expep_sh_chan.srate, 0, [], 0);
%     t.data2=hilbert(squeeze(EEG_expep_sh_ress_f.data));
%     t.data2(:,:,2)=hilbert(squeeze(EEG_expep_sh_chan_f.data));
%     % wavelet from cohen
%     % definte convolution parameters
%     n_wavelet     = EEG_expep_sh_ress.pnts;
%     n_data        = EEG_expep_sh_ress.pnts*EEG_expep_sh_ress.trials;
%     n_convolution = n_wavelet+n_data-1;
%     n_conv_pow2   = pow2(nextpow2(n_convolution));
%     centerfreq = 85/6;
%     time    = -EEG_expep_sh_ress.pnts/EEG_expep_sh_ress.srate/2:1/EEG_expep_sh_ress.srate:EEG_expep_sh_ress.pnts/EEG_expep_sh_ress.srate/2-1/EEG_expep_sh_ress.srate;
%     wavelet = exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*((4/(2*pi*centerfreq))^2)))/centerfreq;
%     eegfft = fft(reshape(EEG_expep_sh_ress.data(1,:,:),1,[]),n_conv_pow2);
%     % convolution
%     eegconv = ifft(fft(wavelet,n_conv_pow2).*eegfft);
%     eegconv = eegconv(1:n_convolution);
%     eegconv = reshape(eegconv(floor((EEG_expep_sh_ress.pnts-1)/2):end-1-ceil((EEG_expep_sh_ress.pnts-1)/2)),EEG_expep_sh_ress.pnts,EEG_expep_sh_ress.trials);
%     t.data3 = eegconv;
%    
%     % plot amplitude time course
%     figure;
%     plot(EEG_expep_sh_ress_gabor.times,abs(mean(t.data(:,:,1),2)));
%     hold on
%     plot(EEG_expep_sh_ress_gabor.times,abs(mean(t.data2(:,:,1),2)));
%     plot(EEG_expep_sh_ress_gabor.times,abs(mean(t.data3(:,:,1),2)));
%     legend({'GABOR';'HILBERT';'WAVELET'})
%     
%     % plot coherence
%     t.itpc_data=squeeze(exp(1i*angle(t.data)));
%     t.itpc_data2=squeeze(exp(1i*angle(t.data2)));
%     t.itpc_data3=squeeze(exp(1i*angle(t.data3)));
%     figure;
%     plot(EEG_expep_sh_ress_gabor.times,abs(mean(t.itpc_data(:,:,1),2)));
%     hold on
%     plot(EEG_expep_sh_ress_gabor.times,abs(mean(t.itpc_data2(:,:,1),2)));
%     plot(EEG_expep_sh_ress_gabor.times,abs(mean(t.itpc_data3(:,:,1),2)));
%     legend({'GABOR';'HILBERT';'WAVELET'})
%     
%     % plot angles
%     figure; 
%     plot(EEG_expep_sh_ress_gabor.times,angle(t.data(:,1,1)))
%     hold on
%     plot(EEG_expep_sh_ress_gabor.times,angle(t.data2(:,1,1)))
%     plot(EEG_expep_sh_ress_gabor.times,angle(t.data3(:,1,1)))
%     legend({'GABOR';'HILBERT';'WAVELET'})
%     
%     figure; 
%     plot(1:EEG_expep_sh_ress_gabor.trials,angle(t.data(1420,:,1)))
%     hold on
%     plot(1:EEG_expep_sh_ress_gabor.trials,angle(t.data2(1420,:,1)))
%     plot(1:EEG_expep_sh_ress_gabor.trials,angle(t.data3(1420,:,1)))
%     legend({'GABOR';'HILBERT';'WAVELET'})
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normans function
%     [ttf.TFR,ttf.timeVec,ttf.freqVec] = mwavelTFA(squeeze(EEG_expep_sh_ress.data),85/6,EEG_expep_sh_ress.srate,85/6*2,3);
%     % calculate average phase value (according to Mike X Cohen)
%     ttf.ITPC_comp = nan(1,numel(EEG_expep_sh_ress.times));
%     for i_tmp = 1:numel(ttf.timeVec)
%         ttf.ITPC_comp(i_tmp)=mean(exp(1i*squeeze(ttf.TFR(1,i_tmp,:))));
%     end
%     ttf.prefAngle_diff = nan(1,numel(EEG_expep_sh_ress.times)-1);
%     for i_tmp=1:numel(EEG_expep_sh_ress.times)-1
%         %angle(ttf.ITPC_comp(i_tmp)-ttf.ITPC_comp(i_tmp+1))
%         ttf.prefAngle_diff(i_tmp)=rad2deg(angdiff2(angle(ttf.ITPC_comp(i_tmp)), angle(ttf.ITPC_comp(i_tmp+1))));
%     end
%     ttf.itpc      = abs(ttf.ITPC_comp);
%     ttf.prefAngle = angle(ttf.ITPC_comp);
%     figure; subplot(3,1,1)
%     plot(EEG_expep_sh_ress.times,ttf.itpc)
%     subplot(3,1,2)
%     plot(EEG_expep_sh_ress.times,ttf.prefAngle)
%     subplot(3,1,3)
%     plot(EEG_expep_sh_ress.times(2:end),ttf.prefAngle_diff)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simple hilbert transform approach
%     EEG_expep_sh_ress_f = pop_eegfiltnew(EEG_expep_sh_ress,(85/6)-0.5,(85/6)+0.5, 32*EEG_expep_sh_ress.srate, 0, [], 0);
%     EEG_expep_sh_chan = pop_select(EEG_expep_sh,'channel',t.chan.StandardLabel);
%     EEG_expep_sh_chan_f = pop_eegfiltnew(EEG_expep_sh_chan,(85/6)-0.5,(85/6)+0.5, 32*EEG_expep_sh_chan.srate, 0, [], 0);
%     t.data=hilbert(squeeze(EEG_expep_sh_ress_f.data));
%     t.data(:,:,2)=hilbert(squeeze(EEG_expep_sh_chan_f.data));
%     
%     ttf.ITPC_comp = nan(2,numel(EEG_expep_sh_ress.times));
%     for i_tmp = 1:numel(EEG_expep_sh_ress.times)
%         ttf.ITPC_comp(1,i_tmp)=mean(exp(1i*(t.data(i_tmp,:,1))));
%         ttf.ITPC_comp(2,i_tmp)=mean(exp(1i*(t.data(i_tmp,:,2))));
%     end
%     ttf.prefAngle_diff = nan(2,numel(EEG_expep_sh_ress.times)-1);
%     for i_tmp=1:numel(EEG_expep_sh_ress.times)-1
%         angle(ttf.ITPC_comp(i_tmp)-ttf.ITPC_comp(i_tmp+1))
%         ttf.prefAngle_diff(1,i_tmp)=rad2deg(angdiff2(angle(ttf.ITPC_comp(1,i_tmp)), angle(ttf.ITPC_comp(1,i_tmp+1))));
%         ttf.prefAngle_diff(2,i_tmp)=rad2deg(angdiff2(angle(ttf.ITPC_comp(2,i_tmp)), angle(ttf.ITPC_comp(2,i_tmp+1))));
%     end
%     ttf.itpc         = abs(ttf.ITPC_comp);
%     ttf.prefAngle    = angle(ttf.ITPC_comp);
%     ttf.ersp_ind     = squeeze(mean(abs(t.data),2))';
%     ttf.ersp_evo     = [abs(hilbert(squeeze(mean(EEG_expep_sh_ress_f.data,3))));...
%         abs(hilbert(squeeze(mean(EEG_expep_sh_chan_f.data,3))))];
%     
%     figure; subplot(2,1,1)
%     plot(EEG_expep_sh_chan_f.times,squeeze(mean(EEG_expep_sh_chan_f.data,3)))
%     hold on
%     plot(EEG_expep_sh_ress.times,ttf.ersp_ind(2,:)); plot(EEG_expep_sh_ress.times,ttf.ersp_evo(2,:)); 
%     title(sprintf('ersp %s',t.chan.StandardLabel{1}))
%     subplot(2,1,2)
%     plot(EEG_expep_sh_chan_f.times,squeeze(mean(EEG_expep_sh_ress_f.data,3)))
%     hold on
%      plot(EEG_expep_sh_ress.times,ttf.ersp_ind(1,:)); plot(EEG_expep_sh_ress.times,ttf.ersp_evo(1,:)); 
%     title(sprintf('ersp RESS'))
%     
%     figure; subplot(3,1,1)
%     plot(EEG_expep_sh_ress.times,ttf.itpc); title('itpc')
%     subplot(3,1,2)
%     plot(EEG_expep_sh_ress.times,ttf.prefAngle); title('angle')
%     subplot(3,1,3)
%     plot(EEG_expep_sh_ress.times(2:end),ttf.prefAngle_diff); title('angle difference')
%     legend({'RESS',t.chan.StandardLabel{1}},'Location','North','Orientation','horizontal')
    
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % coherence measures for alpha
%     ttfa.alpha.freqs = {[10 14]; [8 12]};
%     ttfa.alpha.chans = {{'C3';'CP3'}; {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}};
%     for i_coh = 1:numel(ttfa.alpha.freqs)
%         %EEG_expep_sh_noress_t=pop_select( EEG_expep_n_noress,'channel',ttfa.alpha.chans{i_coh});
%         EEG_expep_sh_noress_t=pop_select( EEG_expep_sh_noress,'channel',ttfa.alpha.chans{i_coh});
%         t.freq_i = dsearchn(F.TFAfreqs',ttfa.alpha.freqs{i_coh}');
%         t.data = nan([size(EEG_expep_sh_noress_t.data), diff(t.freq_i)+1]);
%         for i_freq = 1:diff(t.freq_i)+1
%             [EEG_expep_sh_noress_t_gabor]=eegF_Gabor(EEG_expep_sh_noress_t,F.TFAfreqs(t.freq_i(1)-1+i_freq),0.5,'complex');
%             t.data(:,:,:,i_freq)=EEG_expep_sh_noress_t_gabor.data;
%         end
%         
%         ttf.alpha.itpc_comp{i_coh}  = squeeze(exp(1i*angle(t.data)));
%         ttf.alpha.itpc              = squeeze(abs(mean(ttf.alpha.itpc_comp{i_coh},3)));
%         %     figure; plot(EEG_expep_sh_noress_t_gabor.times,squeeze(ttf.alpha.itpc(2,:,:)))
%         %     figure; plot(EEG_expep_sh_noress_t_gabor.times,squeeze(mean(mean(ttf.alpha.itpc,1),3)))
%         ttf.prefAngle    = angle(ttf.alpha.itpc_comp{i_coh});
%         ttf.alpha.prefAngle_m    = squeeze(angle(mean(mean(mean(ttf.alpha.itpc_comp{i_coh},3),1),4)));
%         %         figure; plot(EEG_expep_sh_noress_t_gabor.times,squeeze(ttf.alpha.prefAngle_m))
%         ttf.ersp_ind     = squeeze(mean(abs(t.data),2));
%         
%     end



    Cohere.itc = ttf.itpc;
    Cohere.itc_comp = ttf.itpc_comp;
    Cohere.itcphase = ttf.prefAngle;
    Cohere.SSVEP_itcphasediff = ttf.prefAngle_diff;
    Cohere.ersp_ind = ttf.ersp_ind;
    Cohere.ersp_evo = ttf.ersp_evo;
    Cohere.freqs = 85/6;
    Cohere.times = EEG_expep_sh_ress.times;
    Cohere.data_index = {'RESS component';sprintf('electrode %s', t.chan.StandardLabel{1})};
    Cohere.filterspecs = 'pop_eegfiltnew(EEG_expep_sh_ress,(85/6)-0.5,(85/6)+0.5, 32*EEG_expep_sh_ress.srate, 0, [], 0);';
    Cohere.method = 'phase coherence based on Gabor transform';
    Cohere.trial_number = EEG_expep_sh_ress.trials;
    Cohere.response_times_n = art.resp_time_n;
    Cohere.response_times_sh = art.resp_time_sh;
    
    
    
    %%%%%%%%%%%%%%%
    % different approach
    % fit sine signal to baseline, extrapolate this and and point by point compare it with
    sinf.sin.f = 85/6;                          % frequency
%     sinf.sin.f = 85/6+7;                          % frequency
    sinf.filt_range = 5;
    sinf.filt_order = 2*EEG_expep_sh_ress.srate; %32*EEG_expep_sh_ress.srate
    EEG_expep_sh_chan = pop_select(EEG_expep_sh,'channel',t.chan.StandardLabel);
    EEG_expep_sh_ress_f = pop_eegfiltnew(EEG_expep_sh_ress,sinf.sin.f-sinf.filt_range,sinf.sin.f+sinf.filt_range, sinf.filt_order, 0, [], 0); %32*EEG_expep_sh_ress.srate
    EEG_expep_sh_chan_f = pop_eegfiltnew(EEG_expep_sh_chan,sinf.sin.f-sinf.filt_range,sinf.sin.f+sinf.filt_range, sinf.filt_order, 0, [], 0);
    EEG_expep_sh_ress_f_m = pop_select(EEG_expep_sh_ress_f,'trial',1); EEG_expep_sh_ress_f_m.data = mean(EEG_expep_sh_ress_f.data,3);
    EEG_expep_sh_chan_f_m = pop_select(EEG_expep_sh_chan_f,'trial',1); EEG_expep_sh_chan_f_m.data = mean(EEG_expep_sh_chan_f.data,3);
    EEG_expep_sh_ress_m = pop_select(EEG_expep_sh_ress,'trial',1); EEG_expep_sh_ress_m.data = mean(EEG_expep_sh_ress.data,3);
    EEG_expep_sh_chan_m = pop_select(EEG_expep_sh_chan,'trial',1); EEG_expep_sh_chan_m.data = mean(EEG_expep_sh_chan.data,3);
    EEG_expep_sh_ress_m_f = pop_eegfiltnew(EEG_expep_sh_ress_m,sinf.sin.f-sinf.filt_range,sinf.sin.f+sinf.filt_range, sinf.filt_order, 0, [], 0); %32*EEG_expep_sh_ress.srate
    EEG_expep_sh_chan_m_f = pop_eegfiltnew(EEG_expep_sh_chan_m,sinf.sin.f-sinf.filt_range,sinf.sin.f+sinf.filt_range, sinf.filt_order, 0, [], 0);
     
%     figure; subplot(2,1,1); plot(EEG_expep_sh_ress.times,[EEG_expep_sh_ress.data(1,:,1);EEG_expep_sh_ress_f.data(1,:,1); EEG_expep_sh_ress_m.data(1,:,1); EEG_expep_sh_ress_m_f.data(1,:,1)])
%     legend({'raw 1trial';'filt 1trial'; 'raw mean'; 'filt mean'}); title('RESS')
%     subplot(2,1,2); plot(EEG_expep_sh_chan.times,[EEG_expep_sh_chan.data(1,:,1);EEG_expep_sh_chan_f.data(1,:,1); EEG_expep_sh_chan_m.data(1,:,1); EEG_expep_sh_chan_m_f.data(1,:,1)])
%     legend({'raw 1trial';'filt 1trial'; 'raw mean'; 'filt mean'}); title(t.chan.StandardLabel{1})
%     
%     figure; subplot(2,1,1); plot(EEG_expep_sh_ress.times,[EEG_expep_sh_ress_f_m.data(1,:,1);EEG_expep_sh_ress_m_f.data(1,:,1); EEG_expep_sh_ress_f_m.data(1,:,1)-EEG_expep_sh_ress_m_f.data(1,:,1)])
%     legend({'filt mean';'mean filt';'diff'}); title('RESS')
%     subplot(2,1,2); plot(EEG_expep_sh_ress.times,[EEG_expep_sh_chan_f_m.data(1,:,1);EEG_expep_sh_chan_m_f.data(1,:,1); EEG_expep_sh_chan_f_m.data(1,:,1)-EEG_expep_sh_chan_m_f.data(1,:,1)])
%     legend({'filt mean';'mean filt';'diff'}); title(t.chan.StandardLabel{1})
    
%     figure; histogram(EEG_expep_sh_chan_f_m.data(1,:,1)-EEG_expep_sh_chan_m_f.data(1,:,1),1000)
    
    % fit sinusoid to baseline window
    sinf.sin.cycle_pnts = ceil(EEG_expep_sh_ress_m_f.srate/sinf.sin.f); % number of points for one cycle
    sinf.baseline = [-3250 -2750];              % baseline window
    sinf.baseline_i = dsearchn(EEG_expep_sh_ress_m_f.times',sinf.baseline');
    sinf.sin.T = 1/EEG_expep_sh_ress_m_f.srate;  % sampling period
%     sinf.sin.L = diff(sinf.baseline_i);         % length of signal
    sinf.sin.L = EEG_expep_sh_ress_m_f.pnts + sinf.sin.cycle_pnts;         % length of signal
    sinf.sin.t =  (0:sinf.sin.L-1)*sinf.sin.T;  % Time vector
    sinf.sin.amp_ress_m = (max(EEG_expep_sh_ress_m_f.data) - min(EEG_expep_sh_ress_m_f.data))/2;
    sinf.sin.amp_chan_m = (max(EEG_expep_sh_chan_m_f.data) - min(EEG_expep_sh_chan_m_f.data))/2;
    sinf.sin.amp_ress_st = (max(squeeze(EEG_expep_sh_ress_f.data)) - min(squeeze(EEG_expep_sh_ress_f.data)))/2;
    sinf.sin.amp_chan_st = (max(squeeze(EEG_expep_sh_chan_f.data)) - min(squeeze(EEG_expep_sh_chan_f.data)))/2;
    sinf.sin.signal_ress_m = sinf.sin.amp_ress_m*sin(2*pi*sinf.sin.f*sinf.sin.t);
    sinf.sin.signal_chan_m = sinf.sin.amp_chan_m*sin(2*pi*sinf.sin.f*sinf.sin.t);
    sinf.sin.signal_ress_st = sinf.sin.amp_ress_st'.*sin(2*pi*sinf.sin.f*sinf.sin.t);
    sinf.sin.signal_chan_st = sinf.sin.amp_chan_st'.*sin(2*pi*sinf.sin.f*sinf.sin.t);
    sinf.bin_width_t = 500;                      % width of bin for extracting lag; in ms
    sinf.bin_width_pnt = ceil(sinf.bin_width_t/1000*EEG_expep_sh_ress_m_f.srate);
    
%     figure; subplot(2,1,1); plot(EEG_expep_sh_ress.times,EEG_expep_sh_ress_m_f.data(:,:,1)); hold on;
%     plot(EEG_expep_sh_ress.times(sinf.baseline_i(1):sinf.baseline_i(2)-1),sinf.sin.signal_ress_m(sinf.baseline_i(1):sinf.baseline_i(2)-1));
%     subplot(2,1,2); plot(EEG_expep_sh_chan.times,EEG_expep_sh_chan_m_f.data(:,:,1)); hold on;
%     plot(EEG_expep_sh_chan.times(sinf.baseline_i(1):sinf.baseline_i(2)-1),sinf.sin.signal_chan_m(sinf.baseline_i(1):sinf.baseline_i(2)-1));
%     
    % align sinusoid to SSVEP signal in baseline time window
    
    
    sinf.sin.bsl_sumofsq_ress_m = []; sinf.sin.bsl_sumofsq_chan_m = [];sinf.sin.bsl_sumofsq_ress_st = []; sinf.sin.bsl_sumofsq_chan_st = [];
    for i_shift = 1:sinf.sin.cycle_pnts
        sinf.sin.signal_ress_m_sh = sinf.sin.amp_ress_m*sin(2*pi*sinf.sin.f*(sinf.sin.t+((i_shift-1)*sinf.sin.T))); % sinusoidal signal RESS
        sinf.sin.bsl_sumofsq_ress_m(i_shift) = sum((sinf.sin.signal_ress_m_sh(sinf.baseline_i(1):sinf.baseline_i(2)-1)-...
            EEG_expep_sh_ress_m_f.data(:,sinf.baseline_i(1):sinf.baseline_i(2)-1,1)).^2); % sum of squares
        
        sinf.sin.signal_chan_m_sh = sinf.sin.amp_chan_m*sin(2*pi*sinf.sin.f*(sinf.sin.t+((i_shift-1)*sinf.sin.T))); % sinusoidal signal chan
        sinf.sin.bsl_sumofsq_chan_m(i_shift) = sum((sinf.sin.signal_chan_m_sh(sinf.baseline_i(1):sinf.baseline_i(2)-1)-...
            EEG_expep_sh_chan_m_f.data(:,sinf.baseline_i(1):sinf.baseline_i(2)-1,1)).^2); % sum of squares
    end
    % figure; subplot(2,1,1); plot((1:sinf.sin.cycle_pnts)-1,sinf.sin.bsl_sumofsq_ress_m);subplot(2,1,2); plot((1:sinf.sin.cycle_pnts)-1,sinf.sin.bsl_sumofsq_chan_m)
    sinf.sin.bsl_shift_ress_m = find(sinf.sin.bsl_sumofsq_ress_m==min(sinf.sin.bsl_sumofsq_ress_m))-1;
    sinf.sin.signal_ress_m_sh = sinf.sin.amp_ress_m*sin(2*pi*sinf.sin.f*(sinf.sin.t+(sinf.sin.bsl_shift_ress_m*sinf.sin.T)));
    sinf.sin.signal_ress_m_sh = sinf.sin.signal_ress_m_sh(1:EEG_expep_sh_ress_m_f.pnts);
    sinf.sin.bsl_shift_chan_m = find(sinf.sin.bsl_sumofsq_chan_m==min(sinf.sin.bsl_sumofsq_chan_m))-1;
    sinf.sin.signal_chan_m_sh = sinf.sin.amp_chan_m*sin(2*pi*sinf.sin.f*(sinf.sin.t+(sinf.sin.bsl_shift_chan_m*sinf.sin.T)));
    sinf.sin.signal_chan_m_sh = sinf.sin.signal_chan_m_sh(1:EEG_expep_sh_chan_m_f.pnts);
    
    % do it for single trial data
    sinf.sin.bsl_sumofsq_ress_st = []; sinf.sin.bsl_sumofsq_chan_st = []; sinf.sin.bsl_shift_ress_st = []; 
    sinf.sin.signal_ress_st_sh = []; sinf.sin.bsl_shift_chan_st = []; sinf.sin.signal_chan_st_sh = [];
    for i_trial = 1:size(sinf.sin.amp_ress_st,2)
        for i_shift = 1:sinf.sin.cycle_pnts
            t.signal_ress_st_sh = sinf.sin.amp_ress_st(i_trial)*sin(2*pi*sinf.sin.f*(sinf.sin.t+((i_shift-1)*sinf.sin.T))); % sinusoidal signal RESS
            sinf.sin.bsl_sumofsq_ress_st(i_shift) = sum((t.signal_ress_st_sh(sinf.baseline_i(1):sinf.baseline_i(2)-1)-...
                EEG_expep_sh_ress_f.data(:,sinf.baseline_i(1):sinf.baseline_i(2)-1,i_trial)).^2); % sum of squares
            
            t.signal_chan_st_sh = sinf.sin.amp_chan_st(i_trial)*sin(2*pi*sinf.sin.f*(sinf.sin.t+((i_shift-1)*sinf.sin.T))); % sinusoidal signal RESS
            sinf.sin.bsl_sumofsq_chan_st(i_shift) = sum((t.signal_chan_st_sh(sinf.baseline_i(1):sinf.baseline_i(2)-1)-...
                EEG_expep_sh_chan_f.data(:,sinf.baseline_i(1):sinf.baseline_i(2)-1,i_trial)).^2); % sum of squares
        end
%         figure; subplot(2,1,1); plot((1:sinf.sin.cycle_pnts)-1,sinf.sin.bsl_sumofsq_ress_st);subplot(2,1,2); plot((1:sinf.sin.cycle_pnts)-1,sinf.sin.bsl_sumofsq_chan_st)
        sinf.sin.bsl_shift_ress_st(i_trial) = find(sinf.sin.bsl_sumofsq_ress_st==min(sinf.sin.bsl_sumofsq_ress_st))-1;
        t.signal_ress_st_sh = sinf.sin.amp_ress_st(i_trial)*sin(2*pi*sinf.sin.f*(sinf.sin.t+(sinf.sin.bsl_shift_ress_st(i_trial)*sinf.sin.T)))';
        sinf.sin.signal_ress_st_sh(1,:,i_trial) = t.signal_ress_st_sh(1:EEG_expep_sh_ress_f.pnts);
        sinf.sin.bsl_shift_chan_st(i_trial) = find(sinf.sin.bsl_sumofsq_chan_st==min(sinf.sin.bsl_sumofsq_chan_st))-1;
        t.signal_chan_st_sh = sinf.sin.amp_chan_st(i_trial)*sin(2*pi*sinf.sin.f*(sinf.sin.t+(sinf.sin.bsl_shift_chan_st(i_trial)*sinf.sin.T)))';
        sinf.sin.signal_chan_st_sh(1,:,i_trial) = t.signal_chan_st_sh(1:EEG_expep_sh_ress_f.pnts);
    end
        
%     figure;  plot(EEG_expep_sh_ress_m_f.times,EEG_expep_sh_ress_m_f.data(:,:,1)); hold on;
%     plot(EEG_expep_sh_ress_m_f.times,sinf.sin.signal_ress_m_sh);
%     xlim(sinf.baseline+[-100 100])
%     figure;  plot(EEG_expep_sh_chan_m_f.times,EEG_expep_sh_chan_m_f.data(:,:,1)); hold on;
%     plot(EEG_expep_sh_chan_m_f.times,sinf.sin.signal_chan_m_sh);
%     xlim(sinf.baseline+[-100 100])
%     figure; subplot(2,1,1);plot(EEG_expep_sh_ress_f.times,squeeze(sinf.sin.signal_ress_st_sh)); xlim(sinf.baseline+[-100 100]);
%     subplot(2,1,2);plot(EEG_expep_sh_chan_f.times,squeeze(sinf.sin.signal_chan_st_sh)); xlim(sinf.baseline+[-100 100])
%     figure; histogram([sinf.sin.bsl_shift_ress_st-sinf.sin.bsl_shift_ress_m]); hold on; histogram(sinf.sin.bsl_shift_chan_st-sinf.sin.bsl_shift_chan_m)

    % calculate lag between baseline reference signal and mean SSVEP
    t.win = -floor(sinf.bin_width_pnt/2):floor(sinf.bin_width_pnt/2);
    sinf.lag_RESS = []; sinf.lag_chan = [];
    for i_time = 1:numel(sinf.sin.signal_ress_m_sh)
        % first extract data for which a sinusoidal signal should be fitted
        t.pntidx = t.win+i_time;
        t.data_ress_signal = [nan(1,sum(t.pntidx<1))...
            EEG_expep_sh_ress_m_f.data(1,t.pntidx>0 & t.pntidx<=numel(sinf.sin.signal_ress_m_sh))...
            nan(1,sum(t.pntidx>numel(sinf.sin.signal_ress_m_sh)))];
        t.data_ress_sin = [nan(1,sum(t.pntidx<1))...
            sinf.sin.signal_ress_m_sh(1,t.pntidx>0 & t.pntidx<=numel(sinf.sin.signal_ress_m_sh))...
            nan(1,sum(t.pntidx>numel(sinf.sin.signal_ress_m_sh)))];
        t.data_chan_signal = [nan(1,sum(t.pntidx<1))...
            EEG_expep_sh_chan_m_f.data(1,t.pntidx>0 & t.pntidx<=numel(sinf.sin.signal_chan_m_sh))...
            nan(1,sum(t.pntidx>numel(sinf.sin.signal_chan_m_sh)))];
        t.data_chan_sin = [nan(1,sum(t.pntidx<1))...
            sinf.sin.signal_chan_m_sh(1,t.pntidx>0 & t.pntidx<=numel(sinf.sin.signal_chan_m_sh))...
            nan(1,sum(t.pntidx>numel(sinf.sin.signal_chan_m_sh)))];
        % hanning
        t.data_ress_signal_H = t.data_ress_signal.*hann(numel(t.data_ress_signal))';
        t.data_ress_sin_H = t.data_ress_sin.*hann(numel(t.data_ress_sin))';
        t.data_chan_signal_H = t.data_chan_signal.*hann(numel(t.data_chan_signal))';
        t.data_chan_sin_H = t.data_chan_sin.*hann(numel(t.data_chan_sin))';
        % figure; plot(t.data_ress_signal); hold on; plot(t.data_ress_sin);
        % loop across shifts
        t.sumofsq_ress_sig = []; t.sumofsq_chan_sig = []; t.sumofsq_ress_sin = []; t.sumofsq_chan_sin = [];
        for i_shift = 1:sinf.sin.cycle_pnts
            % first create sine signal
            t.sin_signal_ress_sh = sinf.sin.amp_ress_m*sin(2*pi*sinf.sin.f*((0:numel(t.win)-1)*sinf.sin.T+((i_shift-1)*sinf.sin.T)));
            t.sin_signal_chan_sh = sinf.sin.amp_chan_m*sin(2*pi*sinf.sin.f*((0:numel(t.win)-1)*sinf.sin.T+((i_shift-1)*sinf.sin.T)));
%             figure; plot(0:numel(t.win)-1,t.sin_signal_ress_sh ); hold on;
%             plot(0:numel(t.win)-1,t.data_ress_signal); plot(0:numel(t.win)-1,t.data_ress_sin);
%             figure; plot(0:numel(t.win)-1,t.sin_signal_chan_sh ); hold on;
%             plot(0:numel(t.win)-1,t.data_chan_signal); plot(0:numel(t.win)-1,t.data_chan_sin);
            
            % apply hanning window
            t.sin_signal_ress_shH = t.sin_signal_ress_sh.*hann(numel(t.sin_signal_ress_sh))';
            t.sin_signal_chan_shH = t.sin_signal_chan_sh.*hann(numel(t.sin_signal_chan_sh))';
            % figure; plot(0:numel(t.win)-1,t.sin_signal_ress_shH ); hold on
            % figure; plot(0:numel(t.win)-1,t.data_ress_signal); hold on; plot(0:numel(t.win)-1,t.data_ress_sin);plot(0:numel(t.win)-1,t.sin_signal_ress_shH )
            % figure; plot(0:numel(t.win)-1,t.data_ress_signal_H); hold on; plot(0:numel(t.win)-1,t.data_ress_sin_H);plot(0:numel(t.win)-1,t.sin_signal_ress_shH )
            
            % calculate square mean differences
%             t.sumofsq_ress_sig(i_shift) = nansum((t.data_ress_signal-t.sin_signal_ress_shH).^2); % sum of squares
%             t.sumofsq_ress_sin(i_shift) = nansum((t.data_ress_sin-t.sin_signal_ress_shH).^2); % sum of squares
%             t.sumofsq_chan_sig(i_shift) = nansum((t.data_chan_signal-t.sin_signal_chan_shH).^2); % sum of squares
%             t.sumofsq_chan_sin(i_shift) = nansum((t.data_chan_sin-t.sin_signal_chan_shH).^2); % sum of squares
            t.sumofsq_ress_sig(i_shift) = nansum((t.data_ress_signal_H-t.sin_signal_ress_shH).^2); % sum of squares
            t.sumofsq_ress_sin(i_shift) = nansum((t.data_ress_sin_H-t.sin_signal_ress_shH).^2); % sum of squares
            t.sumofsq_chan_sig(i_shift) = nansum((t.data_chan_signal_H-t.sin_signal_chan_shH).^2); % sum of squares
            t.sumofsq_chan_sin(i_shift) = nansum((t.data_chan_sin_H-t.sin_signal_chan_shH).^2); % sum of squares
        end
%         figure; subplot(2,1,1); plot((1:sinf.sin.cycle_pnts)-1,t.sumofsq_ress_sig); hold on; plot((1:sinf.sin.cycle_pnts)-1,t.sumofsq_ress_sin);...
%         subplot(2,1,2); plot((1:sinf.sin.cycle_pnts)-1,t.sumofsq_chan_sig); hold on; plot((1:sinf.sin.cycle_pnts)-1,t.sumofsq_chan_sin)
%         
        % save lag between signal and fitted sinusoid (negative: RESS is leading with earlier peak)
        sinf.lag_RESS_m(i_time) = find(t.sumofsq_ress_sig==min(t.sumofsq_ress_sig))-find(t.sumofsq_ress_sin==min(t.sumofsq_ress_sin));
        
        sinf.lag_chan_m(i_time) = find(t.sumofsq_chan_sig==min(t.sumofsq_chan_sig))-find(t.sumofsq_chan_sin==min(t.sumofsq_chan_sin));
    end
%     figure; plot(EEG_expep_sh_ress.times,sinf.lag_RESS_m); hold on;  plot(EEG_expep_sh_ress.times,sinf.lag_chan_m)
    
    % calculate lag between baseline reference signal and single trial SSVEP
    t.win = -floor(sinf.bin_width_pnt/2):floor(sinf.bin_width_pnt/2);
    sinf.lag_RESS_st = nan(size(sinf.sin.signal_ress_st,2),size(sinf.sin.signal_ress_st,1)); sinf.lag_chan_st =  sinf.lag_RESS_st;
    h_wtb = waitbar(0,'1','Name','single trial phase differences...');
    step = 1;
    for i_time = 1:size(sinf.sin.signal_ress_st,2)
        for i_trial = 1:size(sinf.sin.signal_ress_st,1)
            % first extract data for which a sinusoidal signal should be fitted
            t.pntidx = t.win+i_time;
            t.data_ress_signal = [nan(1,sum(t.pntidx<1))...
                EEG_expep_sh_ress_f.data(1,t.pntidx>0 & t.pntidx<=size(sinf.sin.signal_ress_st_sh,2),i_trial)...
                nan(1,sum(t.pntidx>size(sinf.sin.signal_ress_st_sh,2)))];
            t.data_ress_sin = [nan(1,sum(t.pntidx<1))...
                sinf.sin.signal_ress_st_sh(1,t.pntidx>0 & t.pntidx<=size(sinf.sin.signal_ress_st_sh,2),i_trial)...
                nan(1,sum(t.pntidx>size(sinf.sin.signal_ress_st_sh,2)))];
            t.data_chan_signal = [nan(1,sum(t.pntidx<1))...
                EEG_expep_sh_chan_f.data(1,t.pntidx>0 & t.pntidx<=size(sinf.sin.signal_chan_st_sh,2),i_trial)...
                nan(1,sum(t.pntidx>size(sinf.sin.signal_chan_st_sh,2)))];
            t.data_chan_sin = [nan(1,sum(t.pntidx<1))...
                sinf.sin.signal_chan_st_sh(1,t.pntidx>0 & t.pntidx<=size(sinf.sin.signal_chan_st_sh,2),i_trial)...
                nan(1,sum(t.pntidx>size(sinf.sin.signal_chan_st_sh,2)))];
            % figure; plot(t.pntidx, [t.data_ress_signal; t.data_ress_sin])
            % hanning
            t.data_ress_signal_H = t.data_ress_signal.*hann(numel(t.data_ress_signal))';
            t.data_ress_sin_H = t.data_ress_sin.*hann(numel(t.data_ress_sin))';
            t.data_chan_signal_H = t.data_chan_signal.*hann(numel(t.data_chan_signal))';
            t.data_chan_sin_H = t.data_chan_sin.*hann(numel(t.data_chan_sin))';
            % figure; plot(t.data_ress_signal); hold on; plot(t.data_ress_sin);
            % loop across shifts
            t.sumofsq_ress_sig = []; t.sumofsq_chan_sig = []; t.sumofsq_ress_sin = []; t.sumofsq_chan_sin = [];
            for i_shift = 1:sinf.sin.cycle_pnts
                try waitbar(step/(sinf.sin.cycle_pnts*size(sinf.sin.signal_ress_st,1)*size(sinf.sin.signal_ress_st,2)),h_wtb,...
                        sprintf('sub %1.0f of %1.0f | time: %1.0fms',...
                        i_sub, numel(F.Subjects2Use), EEG_expep_sh_ress_f.times(i_time)))
                end
                % first create sine signal
                t.sin_signal_ress_sh = sinf.sin.amp_ress_st(i_trial)*sin(2*pi*sinf.sin.f*((0:numel(t.win)-1)*sinf.sin.T+((i_shift-1)*sinf.sin.T)));
                t.sin_signal_chan_sh = sinf.sin.amp_chan_st(i_trial)*sin(2*pi*sinf.sin.f*((0:numel(t.win)-1)*sinf.sin.T+((i_shift-1)*sinf.sin.T)));
%                 figure; plot(0:numel(t.win)-1,t.sin_signal_ress_sh ); hold on;
%                 plot(0:numel(t.win)-1,t.data_ress_signal); plot(0:numel(t.win)-1,t.data_ress_sin);
%                 figure; plot(0:numel(t.win)-1,t.sin_signal_chan_sh ); hold on;
%                 plot(0:numel(t.win)-1,t.data_chan_signal); plot(0:numel(t.win)-1,t.data_chan_sin);
                
                % apply hanning window
                t.sin_signal_ress_shH = t.sin_signal_ress_sh.*hann(numel(t.sin_signal_ress_sh))';
                t.sin_signal_chan_shH = t.sin_signal_chan_sh.*hann(numel(t.sin_signal_chan_sh))';
%                 figure; plot(0:numel(t.win)-1,t.sin_signal_ress_shH ); hold on
%                 figure; plot(0:numel(t.win)-1,t.data_ress_signal); hold on; plot(0:numel(t.win)-1,t.data_ress_sin);plot(0:numel(t.win)-1,t.sin_signal_ress_shH )
%                 figure; plot(0:numel(t.win)-1,t.data_ress_signal_H); hold on; plot(0:numel(t.win)-1,t.data_ress_sin_H);plot(0:numel(t.win)-1,t.sin_signal_ress_shH )
%                 
                % calculate square mean differences
                %             t.sumofsq_ress_sig(i_shift) = nansum((t.data_ress_signal-t.sin_signal_ress_shH).^2); % sum of squares
                %             t.sumofsq_ress_sin(i_shift) = nansum((t.data_ress_sin-t.sin_signal_ress_shH).^2); % sum of squares
                %             t.sumofsq_chan_sig(i_shift) = nansum((t.data_chan_signal-t.sin_signal_chan_shH).^2); % sum of squares
                %             t.sumofsq_chan_sin(i_shift) = nansum((t.data_chan_sin-t.sin_signal_chan_shH).^2); % sum of squares
                t.sumofsq_ress_sig(i_shift) = nansum((t.data_ress_signal_H-t.sin_signal_ress_shH).^2); % sum of squares
                t.sumofsq_ress_sin(i_shift) = nansum((t.data_ress_sin_H-t.sin_signal_ress_shH).^2); % sum of squares
                t.sumofsq_chan_sig(i_shift) = nansum((t.data_chan_signal_H-t.sin_signal_chan_shH).^2); % sum of squares
                t.sumofsq_chan_sin(i_shift) = nansum((t.data_chan_sin_H-t.sin_signal_chan_shH).^2); % sum of squares
                step = step+1;
            end
            %             figure; subplot(2,1,1); plot((1:sinf.sin.cycle_pnts)-1,t.sumofsq_ress_sig); hold on; plot((1:sinf.sin.cycle_pnts)-1,t.sumofsq_ress_sin);...
            %                 subplot(2,1,2); plot((1:sinf.sin.cycle_pnts)-1,t.sumofsq_chan_sig); hold on; plot((1:sinf.sin.cycle_pnts)-1,t.sumofsq_chan_sin)
            %
            % save lag between signal and fitted sinusoid (negative: RESS is leading with earlier peak)
            sinf.lag_RESS(i_time,i_trial) = find(t.sumofsq_ress_sig==min(t.sumofsq_ress_sig))-find(t.sumofsq_ress_sin==min(t.sumofsq_ress_sin));
            
            sinf.lag_chan(i_time,i_trial) = find(t.sumofsq_chan_sig==min(t.sumofsq_chan_sig))-find(t.sumofsq_chan_sin==min(t.sumofsq_chan_sin));
            figure; plot(sinf.lag_RESS(1:end-1,:))
        end
    end
    figure; plot(EEG_expep_sh_ress.times,sinf.lag_RESS); hold on;  plot(EEG_expep_sh_ress.times,sinf.lag_chan)
    
    
    

    
    %%%%%%%% save data
    ttf.itpc_comp=squeeze(exp(1i*angle(t.data)));
        
    ttf.itpc         = squeeze(abs(mean(ttf.itpc_comp,2)));
    ttf.prefAngle    = angle(ttf.itpc_comp);
    ttf.ersp_ind     = squeeze(mean(abs(t.data),2));
    ttf.ersp_evo     = [abs(EEG_expep_sh_ress_m_gabor.data(:,:,1))',abs(EEG_expep_sh_chan_m_gabor.data(:,:,1))'];
    
    ttf.prefAngle_diff = nan(numel(EEG_expep_sh_ress.times)-1,EEG_expep_sh_ress_gabor.trials,2);
    for i_tmp=1:numel(EEG_expep_sh_ress.times)-1
        ttf.prefAngle_diff(i_tmp,:,:)=rad2deg(angdiff2(angle(ttf.itpc_comp(i_tmp,:,:)), angle(ttf.itpc_comp(i_tmp+1,:,:))));
    end
    
    
    
    %% save
    Cohere.electrodes = EEG_expep_sh.chanlocs;
    Cohere.electrodes_RESS = EEG_expep_sh_ress.chanlocs;
    Cohere.params.srate = EEG_expep_n.srate;
    Cohere.savetime=datestr(now);
    Cohere.blink_free_win = art.win_blinks;
    Cohere.Gabor_FWHM = F.TFA_FWHM;
    fprintf(1,'|| saving file ||  %s\\VP%02.0f_exp_cohere.mat ||\n', F.PathOut,F.Subjects2Use(i_sub)')
    save(sprintf('%s\\VP%02.0f_exp_cohere.mat', F.PathOut,F.Subjects2Use(i_sub)), 'Cohere')   
    
end
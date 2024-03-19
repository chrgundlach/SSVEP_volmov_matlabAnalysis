%% script in order to calculate regression based on SSVEPs
% based on RESS components



clearvars

%% parameters
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\RAW';
F.PathOut               = 'D:\work\data\SSVEP_volmov\EEG\TFA_lagged_regression';
F.PathInEOG             = 'D:\work\data\SSVEP_volmov\EEG\EOG_Data\';
F.PathInRESS            = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD2';
F.PathInRESS            = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD';
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
% F.Subjects2Use          = [12];
F.Subjects2Use          = [1];
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.VEOGChans             = [];
F.HEOGChans             = [];
F.ChanlocFile           = 'C:\Users\psy05cvd\Dropbox\work\matlab\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.RespTrigger           = {'30'};
F.RespEpoch             = [-5 5]; % in s
F.SSVEPTrigger          = {'15'};
F.BlinkTrigger          = {'50'};
F.BlinkWindow           = [-1000 1000]; % time window in ms around blink for which data will be discarded

F.ExpEpochTrigger       = {'11','12'};
F.eptime                = 359.5; % time of each block in s
F.Blocks                = 6; % number of blocks

F.FrameRate             = 85;
F.TriggerRate           = 6; % SSVEP trigger frequency (e.g. 6 = every 6s)
F.flickframes           = [1 1 1 0 0 0]; % on-off frames for SSVEP

%F.TFAfreqs              = [5:(1/6):40];
F.TFAfreqs              = [4-(1/3):0.5:45];
F.TFAfreqs              = [2-(1/3):0.25:31];

F.CSD_flag              = 1; % 0 = no; 1 = yes


Regress.timewin         = [-2000 2000];
Regress.timeexclude     = [1000]; % esclude first and last 1000ms due to Gabor smearing

%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp.bdf ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
%     EEG=pop_biosig(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)));
    EEG = pop_readbdf(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)), [] ,F.EMGChans(2)+1,[]);
    % pop_eegplot(EEG,1,1,1);
    
    % troubleshooting VP24
%     EEG = pop_loadset('filename','VP24_exp.set','filepath','F:\\work\\data\\SSVEP_volmov\\EEG\\RAW\\');
    
    % load RESS data
    TFA = open(sprintf('%s\\VP%02.0f_exp_tfa.mat',F.PathInRESS,F.Subjects2Use(i_sub)));
    
    % chanlocs
    EEG.chanlocs=pop_chanedit(EEG.chanlocs,'load',{F.ChanlocFile,'filetype','besa (elp)'}); % Load Channel Locations
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    
    % rereference
    EEG = pop_reref( EEG, [],'exclude',[F.EEGChans(end)+1:EEG.nbchan] ); % average
    % pop_eegplot(EEG,1,1,1,1)
    
    
    %% find blinks
    % bipolarize data
    EEG_EOG=eegF_Bipolarize(EEG,{EEG.chanlocs(65:68).labels});
    EEG_EOG=pop_select(EEG_EOG,'channel',{'VEOG','HEOG'});
    EEG_EOG_f = pop_eegfiltnew(EEG_EOG,0.5, 0, 8*EEG_EOG.srate, 0, [], 0);
    %figure; plot(EEG_exp_EOG_f.times,EEG_exp_EOG_f.data(1,:,1));
    [blinks.peaks blinks.peaklocs] = findpeaks(EEG_EOG_f.data(1,:),'MinPeakDistance',EEG_EOG_f.srate/2,'MinPeakHeight',100);
    %figure; plot(EEG_exp_EOG_f.times,EEG_exp_EOG_f.data(1,:,i_bl));
    %hold on; plot(EEG_exp_EOG_f.times(blinks.peaklocs),EEG_exp_EOG_f.data(1,blinks.peaklocs,i_bl),'or')
    
    % write events of blinks
    for i_ev = 1:numel(blinks.peaklocs)
        EEG = pop_editeventvals(EEG,'append',{1 [] [] []},...
            'changefield',{2 'latency' blinks.peaklocs(i_ev)/EEG_EOG_f.srate},... % in seconds
            'changefield',{2 'type' str2num(F.BlinkTrigger{1})});
    end
   
    EEG=eeg_checkset(EEG);
    
    %% use RESS component
    % use spatial filters from RESS data for data without epochs
    EEG_raw = pop_select(EEG,'channel',1:64);
    EEG_raw = pop_resample(EEG_raw, 256);
    EEG_ress = EEG_raw;
    EEG_noress = EEG_raw;
    
    % use spatial filters ofshifted data for not shifted data
    EEG_ress = pop_select( EEG_ress,'channel',1);
    EEG_ress=pop_chanedit(EEG_ress, 'changefield',{1 'labels' 'RESS'},...
        'changefield',{1 'theta' ''},'changefield',{1 'radius' ''},'changefield',...
        {1 'X' ''},'changefield',{1 'Y' ''},'changefield',{1 'Z' ''},'changefield',{1 'sph_theta' ''},...
        'changefield',{1 'sph_phi' ''},'changefield',{1 'sph_radius' ''});
    
    % ress data
    ress_ts = zeros(size(EEG_ress.data));
    for i_tr=1:EEG_ress.trials
        ress_ts(1,:,i_tr) = TFA.TFA.RESS.eigenvectors(:,TFA.TFA.RESS.Comp_selected)'*squeeze(EEG_raw.data(:,:,i_tr));
    end
    EEG_ress.data=ress_ts;
    %figure; plot(EEG_ress.times,EEG_ress.data);
    %pop_eegplot(EEG_ress,1,1,1)
    
    % noress data
    noress_ts = zeros(size(EEG_noress.data));
    for i_tr=1:EEG_noress.trials
        noress_ts(:,:,i_tr)  = TFA.TFA.RESS.ComponentMaps(:,1:end-1) * TFA.TFA.RESS.eigenvectors(:,1:end-1)' * EEG_raw.data(:,:,i_tr);
    end
    EEG_noress.data=noress_ts;
    %pop_eegplot(EEG_noress,1,1,1)
    
    if F.CSD_flag == 1 & i_sub == 1 % calculate CSD matrix
        CSD.chanmat=ExtractMontage('C:\Users\psy05cvd\Dropbox\work\matlab\software\toolboxes\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG_raw.chanlocs.labels}');
        [CSD.G,CSD.H] = GetGH(CSD.chanmat);
    end
    
    
    
    %% epoch data and apply CSD transform
    %     t.t1 = cell2mat({EEG_ress.event(cell2mat({EEG_ress.event.type}) == 11).latency});
    %     t.t2 = cell2mat({EEG_ress.event(cell2mat({EEG_ress.event.type}) == 253).latency});
    %     min(abs(t.t1(1:end-1)-t.t2(2:end)))/EEG_ress.srate;
    
    
    EEG_ress_ep = pop_epoch(EEG_ress, F.ExpEpochTrigger(1), [0 F.eptime], 'epochinfo', 'yes');
    if EEG_ress_ep.trials < 6
        EEG_ress_ep = pop_epoch(EEG_ress, F.ExpEpochTrigger(1), [0 F.eptime-0.5], 'epochinfo', 'yes');
    end
    % pop_eegplot(EEG_ress_ep,1,1,1,1)
    EEG_noress_ep = pop_epoch(EEG_noress, F.ExpEpochTrigger(1), [0 F.eptime], 'epochinfo', 'yes');
    if EEG_noress_ep.trials < 6
        EEG_noress_ep = pop_epoch(EEG_noress, F.ExpEpochTrigger(1), [0 F.eptime-0.5], 'epochinfo', 'yes');
    end
    % pop_eegplot(EEG_noress_ep,1,1,1,1)
    
    
    fprintf(1,'\n###\ncalculating CSD transform\n###\n')
    if F.CSD_flag == 1
        for i_tr = 1:EEG_noress_ep.trials
            % csd of raw data
            EEG_noress_ep.data(:,:,i_tr)= CSDTransform(EEG_noress_ep.data(:,:,i_tr), CSD.G, CSD.H);
        end
    end
    
    %% start of actual regression procedure
    % extract signal to function as a regressor, e.g. SSVEP time course
    % use gabor
    TFA.type = 'gabor';
    TFA.f = F.TFAfreqs;
    TFA.params.gabor_FWHM_freq = 0.5;
    TFA.params.gabor_FWHM_freq = 1;
    TFA.resize = 8; % resize factor
    
    % define regressor signal
    TFA.regressor_freq = F.FrameRate/numel(F.flickframes);
    t.data = nan([floor(EEG_ress_ep.pnts/TFA.resize) EEG_ress_ep.trials  numel(TFA.regressor_freq)]);
    for i_freq = 1:numel(TFA.regressor_freq)
        fprintf(1,'regressor signal:...\n')
        EEG_ress_ep_gab = eegF_Gabor(EEG_ress_ep,TFA.regressor_freq(i_freq), TFA.params.gabor_FWHM_freq);
        %figure; plot(EEG_ress_ep_gab.times,squeeze(EEG_ress_ep_gab.data))
        
        %%%%%%%%%% remove data around blink
        for i_tr = 1:EEG_ress_ep_gab.trials
            t.blinksamples = cell2mat(EEG_ress_ep_gab.epoch(i_tr).eventlatency(...
                cell2mat([EEG_ress_ep_gab.epoch(i_tr).eventtype])==str2num(F.BlinkTrigger{1})))/1000*EEG_ress_ep_gab.srate;
            % loop across blinks to replace data by nan
            t.blinkwindow_i = F.BlinkWindow./1000*EEG_ress_ep_gab.srate;
            for i_blink = 1:numel(t.blinksamples)
                % index values closest to blink.window (even if blink is close to start end)
                t.blinkwindow_i2 = dsearchn((1:EEG_ress_ep_gab.pnts)', (t.blinksamples(i_blink)+t.blinkwindow_i(1))):...
                    dsearchn((1:EEG_ress_ep_gab.pnts)', (t.blinksamples(i_blink)+t.blinkwindow_i(2)));
                % replace values around blink with nans
                EEG_ress_ep_gab.data(1,t.blinkwindow_i2,i_tr) = nan(size(t.blinkwindow_i2));
            end
        end
        % pop_eegplot(EEG_ress_ep_gab,1,1,1,1)
        %%%%%%%%%%%
        
        %%%%%%%%%%%
        % get rid of high and low values
        % figure; histogram(EEG_ress_ep_gab.data(:))
        EEG_ress_ep_gab.data(...
            EEG_ress_ep_gab.data(:) <= (nanmean(EEG_ress_ep_gab.data(:))-(3*nanstd(EEG_ress_ep_gab.data(:)))) |...
            EEG_ress_ep_gab.data(:) >= (nanmean(EEG_ress_ep_gab.data(:))+(3*nanstd(EEG_ress_ep_gab.data(:))))...
            ) = nan;
        %%%%%%%%%%%
        
        t.data(:,:,i_freq) = imresize(squeeze(mean(EEG_ress_ep_gab.data,1)), ...
            [floor(EEG_ress_ep_gab.pnts/TFA.resize) EEG_ress_ep_gab.trials]);
        
        %%% for checking
        %EEG_noress_ep_gab = eegF_Gabor(EEG_noress_ep,TFA.regressor_freq(i_freq), TFA.params.gabor_FWHM_freq);
        %t.data(:,:,i_freq) = imresize(squeeze(mean(EEG_noress_ep_gab.data(29,:,:),1)), ...
        % [floor(EEG_noress_ep_gab.pnts/TFA.resize) EEG_noress_ep_gab.trials]);
        
    end
    TFA.regressor_data = mean(t.data,3);
    TFA.regressor_time = imresize(EEG_ress_ep_gab.times,[1 floor(EEG_ress_ep_gab.pnts/TFA.resize)]);
    % figure; plot(TFA.regressor_time,TFA.regressor_data)
    
    % define source signal and do regression for each frequency (due to memory problems)
    TFA.source_freq = TFA.f;
    
    %% loop across each frequency
    h_wtb = waitbar(0,'1','Name','SSVEP regression...');
    step = 1;
    tic
    fprintf(1,'predictor signal:...\n')
    for i_freq = 1:numel(TFA.source_freq)      
        %% loop across frequencies
        % calculate gabor based TFA-estimation for frequency of interest
        EEG_noress_ep_gab = eegF_Gabor(EEG_noress_ep,TFA.source_freq(i_freq), TFA.params.gabor_FWHM_freq);
%         EEG_noress_ep_gab = eegF_Detrend(EEG_noress_ep_gab);
        %figure; plot(EEG_noress_ep_gab.times,squeeze(EEG_noress_ep_gab.data(29,:,:)))
        % t.tdata = EEG_noress_ep_gab.data(29,:,:); figure; plot(t.tdata(:))
        
        % resize data
        t.source_data = nan([EEG_noress_ep_gab.nbchan floor(EEG_noress_ep_gab.pnts/TFA.resize) EEG_noress_ep_gab.trials]);
        for i_tr = 1:EEG_noress_ep_gab.trials
            t.source_data(:,:,i_tr) = imresize(EEG_noress_ep_gab.data(:,:,i_tr), ...
                [EEG_noress_ep_gab.nbchan floor(EEG_noress_ep_gab.pnts/TFA.resize)]);
        end
        
        %%%%%%%%%%% replace data around blink
        for i_el = 1:size(t.source_data,1)
            t.source_data(i_el,isnan(TFA.regressor_data))=nan(1,sum(sum(isnan(TFA.regressor_data))));
        end
        % figure; plot(t.source_data(1,:,2)); hold on; plot(TFA.regressor_data(:,2).*200)
        %%%%%%%%%%%
        
        % actual regression
        if i_freq == 1
            % define time of source data
            TFA.source_time = imresize(EEG_noress_ep_gab.times,[ 1 floor(EEG_noress_ep_gab.pnts/TFA.resize)]);
            % define time values for which regrssion is done (in ms)
            REG.time = (Regress.timewin(1):1000/(EEG_noress_ep_gab.srate/TFA.resize):Regress.timewin(2));
            % get index of datapoints representing each lag
            REG.pntindex = REG.time/(1000/(EEG_noress_ep_gab.srate/TFA.resize));
            % preallocate betas, p-values (whole model), t-statvalues,R-squared-adjusted
            REG.data_b = nan(EEG_noress_ep_gab.nbchan,numel(TFA.source_freq),numel( REG.time));
            REG.data_p = nan(EEG_noress_ep_gab.nbchan,numel(TFA.source_freq),numel( REG.time));
            REG.data_tstat = nan(EEG_noress_ep_gab.nbchan,numel(TFA.source_freq),numel( REG.time));
            REG.data_Rsqr_adj = nan(EEG_noress_ep_gab.nbchan,numel(TFA.source_freq),numel( REG.time));
            % frequencies
            REG.freq = TFA.source_freq;
            % get channels
            REG.chans = {EEG_noress_ep_gab.chanlocs.labels};
        end
        
        % loop for each lag
        for i_lag = 1:numel(REG.pntindex)
            %% loop across lags
            try waitbar(step/(numel(REG.data_b)),h_wtb,...
                    sprintf('sub %1.0f of %1.0f | freq: %1.2fHz | lag: %1.0fms',...
                    i_sub, numel(F.Subjects2Use), REG.freq(i_freq),REG.time(i_lag)))
            end
            
            % index data
            t.ind_regressor = find(TFA.regressor_time >= (-Regress.timewin(1) + Regress.timeexclude) & ...
                TFA.regressor_time <= TFA.regressor_time(end)-(-Regress.timewin(1) + Regress.timeexclude));
            t.ind_source = t.ind_regressor+REG.pntindex(i_lag);
            
            % extract data across trials
            t.data_regressor = reshape(TFA.regressor_data(t.ind_regressor,:),[],1);
            t.data_source = reshape(t.source_data(:,t.ind_source,:),size(t.source_data,1),[]);
                    
            % graphical checking
            % figure; plot(t.data_source(29,:)); figure; plot(t.data_regressor)
            %figure; scatter(t.data_regressor,t.data_source(29,:),'.')
            %figure; scatter(zscore(t.data_regressor),zscore(t.data_source(29,:)),'.')
            
            % alternative approach
            % regress [,...stats] = regress(data to be predicted, predictor);
            % stats = R2 statistic, F statistic, p value,  estimate of the error variance
            %[regr.b,regr.bint,regr.r,regr.rint,regr.stats] = ...
            %    regress(t.data_regressor,[ones(numel(t.data_regressor),1) t.data_source(29,:)']);
            
            % loop across electrodes
            for i_el = 1:size(t.data_source,1)
                %% loop across electrodes
                % index nans in any of the two
                t.nanindex = ~(isnan(t.data_source(i_el,:)') | isnan(t.data_regressor));
                
                % linear regeression
                regr.lm = fitlm(t.data_source(i_el,t.nanindex ),t.data_regressor(t.nanindex ),'linear');
                %figure; scatter(t.data_regressor(t.nanindex ),t.data_source(i_el,t.nanindex ),'.')
                %figure; scatter(zscore(t.data_regressor(t.nanindex )),zscore(t.data_source(i_el,t.nanindex )),'.')
                [rho,pval]=corr(t.data_source(i_el,t.nanindex )',t.data_regressor(t.nanindex ));
                
                % get estimate for electrode, frequence, lag
                REG.data_b(i_el,i_freq,i_lag) = regr.lm.Coefficients.Estimate(2);
                
                % get tstat for electrode, frequence, lag
                REG.data_tstat(i_el,i_freq,i_lag) = regr.lm.Coefficients.tStat(2);
                
                % get p-value for electrode, frequence, lag
                REG.data_p(i_el,i_freq,i_lag) = regr.lm.Coefficients.pValue(2);
                
                % get R-squared for electrode, frequence, lag
                REG.data_Rsqr_adj(i_el,i_freq,i_lag) = regr.lm.Rsquared.Adjusted;
                
                % first graphical checking
%                 figure; 
%                 subplot(3,1,1); plot( log10(REG.data_p(:,i_freq,i_lag))); hline(log10(.05),'r'); title('p-vals')
%                 subplot(3,1,2); plot( REG.data_Rsqr_adj(:,i_freq,i_lag)); title('R^2')
%                 subplot(3,1,3); plot( REG.data_tstat(:,i_freq,i_lag)); title('t-stat')
                
                step = step + 1;
            end            
        end
%         % graphical checking for different lags and electrodes
%         figure; plot(REG.time,squeeze(log10(REG.data_p(:,i_freq,:)))'); hline(log10(.05),'r'); title('p-vals')
%         figure; plot(REG.time,squeeze(REG.data_Rsqr_adj(:,i_freq,:))');  title('R^2')
%         figure; plot(REG.time,squeeze(REG.data_tstat(:,i_freq,:))');  title('t-stat')
                
    end
    close(h_wtb)
    
    % save relevant data
    REG.chanlocs = EEG_noress_ep_gab.chanlocs;
    REG.params.savetime = datestr(now);
    REG.params.elapsedtime = toc;
    REG.params.TFA_type = TFA.type;
    REG.params.TFA_gabor_FWHM_freq = TFA.params.gabor_FWHM_freq;
    REG.params.TFA_resize = TFA.resize;
    REG.params.TFA_regressor_freq = TFA.regressor_freq;
    REG.params.TFA_predictor_freq = TFA.source_freq;
    REG.params.RESS = TFA.TFA.RESS;
    
    fprintf(1,'|| saving file ||  %s\\VP%02.0f_exp_lagreg_noblinks.mat ||\n', F.PathOut,F.Subjects2Use(i_sub)')
    if ~exist(F.PathOut); mkdir(F.PathOut); end
    save(sprintf('%s\\VP%02.0f_exp_lagreg_withSSVEP_noblinks.mat', F.PathOut,F.Subjects2Use(i_sub)), 'REG')   
    
    %% graphical checking
%      % version 1
%     %%% start
%     pl.elec2plot = {'O1';'Oz';'O2';'POz'}; pl.ROI = 'visual central';
% %     pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.ROI = 'vis alpha/beta';
% %     pl.elec2plot = {'C3';'CP3'}; pl.ROI = 'mot alpha';
% %     pl.elec2plot = {'FP1';'FPz';'FP2'}; pl.ROI = 'frontal eyes';
% %     pl.elec2plot = {'F1';'Fz';'F2';'FC1';'FCz';'FC2'}; pl.ROI = 'frontal';
%     pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({REG.chanlocs.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
%     
%     pl.plotdata = []; pl.clims = []; pl.titadd = {};
%     pl.plotdata(:,:,1) = squeeze(mean(REG.data_b(pl.elec2plot_i,:,:),1)); pl.clims(1,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,1)))); pl.titadd{1} ='betas';
%     pl.plotdata(:,:,2) = squeeze(mean(REG.data_tstat(pl.elec2plot_i,:,:),1)); pl.clims(2,:) = [-1 1]*max(max(abs(pl.plotdata(:,:,2)))); pl.titadd{2} ='t-stat';
%     pl.plotdata(:,:,3) = squeeze(abs(log10(mean(REG.data_p(pl.elec2plot_i,:,:),1)))); pl.titadd{3} ='log p-val';
%     pl.tdata = pl.plotdata(:,:,3); pl.clims(3,:) = [0 1]*max(max(abs(pl.tdata(~isinf(pl.tdata))))); 
%     pl.plotdata(:,:,4) = squeeze(mean(REG.data_Rsqr_adj(pl.elec2plot_i,:,:),1)); pl.clims(4,:) = [0 1]*max(max(abs(pl.plotdata(:,:,4)))); pl.titadd{4} ='r^2 adj';
%     
%     figure;
%     for i_d = 1:4
%         subplot(2,2,i_d)
%         imagesc(REG.time,REG.freq,pl.plotdata(:,:,i_d),pl.clims(i_d,:))
%         set(gca,'YDir','normal')
%         title(sprintf('time lagged regression on SSVEP signal | %s\n%s',  pl.ROI,  pl.titadd{i_d}), 'FontSize',8)
%         colorbar
%         xlabel('lag in ms')
%         ylabel('frequency in Hz')
%         hline(14.16667,'m')
%         set(gca,'FontSize',8)
%         colormap('jet')
%     end
    %%% end
    
    %% version 2 (in topography)
%     pl.plotdata = []; pl.clims = []; pl.titadd = {};
%     pl.plotdata(:,:,:,1) = REG.data_b; pl.clims(1,:) = [-1 1]*max(max(max(abs(pl.plotdata(:,:,:,1))))); pl.titadd{1} ='betas';
%     pl.plotdata(:,:,:,2) = REG.data_tstat; pl.clims(2,:) = [-1 1]*max(max(max(abs(pl.plotdata(:,:,:,2))))); pl.titadd{2} ='t-stat';
%     pl.plotdata(:,:,:,3) = abs(log10(REG.data_p)); pl.titadd{3} ='log p-val';
%     pl.tdata = pl.plotdata(:,:,:,3); pl.clims(3,:) = [0 150]; %pl.clims(3,:) = [0 1]*max(max(max(abs(pl.tdata(~isinf(pl.tdata)))))); 
%     pl.plotdata(:,:,:,4) = REG.data_Rsqr_adj; pl.titadd{4} ='r^2 adj';
%     pl.clims(4,:) = [0 0.01]; %pl.clims(4,:) = [0 1]*max(max(max(abs(pl.plotdata(:,:,:,4))))); 
%      
%     for i_d = 1:4
%         plot_data_topoarray(REG.chanlocs, pl.plotdata(:,:,:,i_d),'TFA','times',REG.time,...
%             'freqs',REG.freq,'title', pl.titadd{i_d},'clim', pl.clims(i_d,:),'ylim_f',[REG.freq(1) 25])
%     end
    


    
    
    
 
end




function [RESS] = RESS_Calculate_Component(EEG_in, freqs, varargin)
%   RESS_CALCULATE_COMPONENT(EEG_in, freqs, varargin)
%   Function to calculate RESS component
%   based on publication
%    "Rhythmic entrainment source separation:
%     Optimizing analyses of neural responses to rhythmic sensory stimulation"
%   from Mike X Cohen and Rasa Gulbinaite
%
%   Cohen, M. X., & Gulbinaite, R. (2016). Rhythmic entrainment source separation: 
%   Optimizing analyses of neural responses to rhythmic sensory stimulation. NeuroImage, 
%   147, 43–56. https://doi.org/10.1016/j.neuroimage.2016.11.036
%
% requires:
%       EEG_in based on eeglab format
%       freqs
%
% possible input
%       'trigger', {'S  1';'S 12'};         uses only specified trigger for component calculation
%       'time', [500 1500];                 uses only specified time for component calculations
%       'peakwidth', [.5];                  FWHM at peak frequency
%       'neighfreqdist', [1];               distance of neighboring frequencies away from peak frequency +/- in Hz;
%       'neighfreqwidth', [1];              FWHM of the neighboring frequencies
%       'neighparams', {[-1 +1; 1 1];[-2 +1; 1 1]}
%                                           manually specify neighbor/reference parameters separately for freqs and low and high reference signals
%                                           {[l_dist r_dist; l_width r_width]; [l_dist r_dist; l_width r_width]}
%       'diagnostics', 'on';                plots graphical diagnostics
%       'broadband', true                   use broadband signal as noise matrix
%
% output
%       'EEG_out'                           contains only RESS component
%       'EEG_out_noRess'                    contains all channels, with RESS component 'subtracted'
%       'EEG_out_withRESS'                  contains RESS component as 'additional channel'
%       'Mean'                              averaged data



% to be implemented
%   - SNR estimation for each component                            done!
%   - R spectrum paramters differently shifted?

%C.Gundlach 2016 (chr.gundlach@gmail.com)
% edited:   2017 | added different output format, C.Gundlach, v.2017-12-19
%           2018 | added possibility to differentially specify reference
%           signals for high and low v.2018-06-05
%           2018 | use of broadband-signal for noise matrix implemented v.2018-06-19
%           2018 | minor bugs fixed (i.e. components for all trials, independent of trigger selection) v.2018-06-22
%           2018 | major bug fixed (width, frequency selection was messed) v.2018-09-12
%           2020 | added regularization of noise for increased fidelity of estimation v.2020-06-18
%                   see Gulbinaite et al. (2019, NeuroImage) https://doi.org/10.1016/j.neuroimage.2019.116146


% current version: v.2020-06-18

%% presets
if nargin <=2
    help RESS_Calculate_Component;
    return
end

trials2use = true(1,EEG_in.trials);         % index of trials used for component calculation
time2use = true(1,EEG_in.pnts);             % index of timepoints to be used for calculation
peakwidt  = .5;                             % FWHM at peak frequency
neighfreq = 1;                              % distance of neighboring frequencies away from peak frequency, +/- in Hz
neighwidt = 1;                              % FWHM of the neighboring frequencies
broadband = false;                          % use broadband signal instead of filter?
epochvect = ones(1,EEG_in.trials);          % epoch vector
diagnosticsflag = 0;                        % by default don't output diagnostics
timelims = [EEG_in.xmin EEG_in.xmax]*1000; % default time limits dependent on input data
neighparams = repmat({[-neighfreq neighfreq; neighwidt neighwidt]},numel(freqs),1);
regparam = .01;                             % noise cov regularization suggested: .01
ress_v = 'v.2020-06-18';                    % current version



%% check for input
if ~(round(numel(varargin)/2) == numel(varargin)/2)
    error('Odd number of input arguments??')
end

for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    if ~isstr(Param)
        error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
        case 'trigger' % define which triggers 2 use
            epochvect = zeros(1,EEG_in.trials);
            for i=1:EEG_in.trials
                [~,t0] = min(abs(cell2mat(EEG_in.epoch(i).eventlatency)));
                try epochvect(i) = EEG_in.epoch(i).eventtype{t0};
                catch
                    epochvect(i) = str2num(EEG_in.epoch(i).eventtype{t0});
                end
            end
            trials2use=ismember(epochvect,cell2mat(Value));
        case 'time' % define time window
            t.index = dsearchn(EEG_in.times',Value');
            time2use = false(1,EEG_in.pnts);
            time2use(t.index(1):t.index(2))=true;
            timelims = Value;
        case 'peakwidth' % define peak freq width
            peakwidt=Value;
        case 'neighfreqdist' % define distance of neighboring frequencies
            neighfreq=Value;
            neighparams = repmat({[-neighfreq neighfreq;  neighwidt neighwidt]},numel(freqs),1);
        case 'neighfreqwidth' % define width of neighboring frequencies
            neighwidt=Value;
            neighparams = repmat({[-neighfreq neighfreq; neighwidt neighwidt]},numel(freqs),1);
        case 'neighparams'
            neighparams=Value;
        case 'broadband' % maximization of signal agains broadband
            broadband = true;
            if Value
                neighparams = 'broadband';
            end
        case 'diagnostics' % set diagnostics
            if strcmpi(Value,'on')
                diagnosticsflag = 1;
            end
    end
end

%% start calculation
% loop for frequencies
fprintf('\ncalculating RESS-component(s) for %1.0f frequencies...\n',numel(freqs))
for i_freq = 1:numel(freqs)
    % test
    %     EEG_in = pop_resample(EEG_in, 256);
    
    
    % FFT parameters
    nfft = ceil(EEG_in.srate/.1); % .1 Hz resolution
    
    % extract EEG data
    data  = EEG_in.data(:,:,trials2use);
    dataX = mean(abs(fft(data(:,time2use,:),nfft,2)/sum(time2use)).^2,3);
    dataXm = abs(fft(mean(data(:,time2use,:),3),nfft,2)/sum(time2use)).^2;
    hz    = linspace(0,EEG_in.srate,nfft);
    % figure; plot(hz,mean(dataX)); xlim([0 60]);
    % figure; plot(hz,mean(dataXm)); xlim([0 60]);
    % figure; plot(hz,dataX(29,:)); xlim([0 60]);
    % figure; plot(hz,dataXm(29,:)); xlim([0 60]);
    
    % topoplot of
    % tfix = dsearchn(hz',freqs(i_freq));
    % figure; topoplot(dataX(:,tfix), EEG_in.chanlocs, 'maplimits',[0 max(dataX(:,tfix))]); colormap('hot'); colorbar;
    % figure; topoplot(dataXm(:,tfix), EEG_in.chanlocs, 'maplimits',[0 max(dataXm(:,tfix))]); colormap('hot'); colorbar;
    
    
    % compute covariance matrix at peak frequency
    fdatAt = filterFGx(data, EEG_in.srate,freqs(i_freq),peakwidt);
    %     figure; plot(EEG_in.times,mean(mean(data(29,:,:),3),1)); title('signal')
    %     figure; plot(EEG_in.times,mean(mean(fdatAt(29,:,:),3),1)); title('signal filtered')
    fdatAt = reshape( fdatAt(:,time2use,:), EEG_in.nbchan,[] );
    fdatAt = bsxfun(@minus,fdatAt,mean(fdatAt,2));
    covAt  = (fdatAt*fdatAt')/sum(time2use);
    % figure; imagesc(covAt); axis square
    
    % compute covariance matrix for lower neighbor
    if ~broadband
        fdatLo = filterFGx(data,EEG_in.srate,freqs(i_freq)+neighparams{i_freq}(1,1),neighparams{i_freq}(2,1));
    else
        fdatLo = data;
    end
    %     figure; plot(EEG_in.times,mean(mean(fdatLo(29,:,:),3),1)); title('noise 1')
    fdatLo = reshape( fdatLo(:,time2use,:), EEG_in.nbchan,[] );
    fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
    covLo  = (fdatLo*fdatLo')/sum(time2use);
    % figure; imagesc(covLo); axis square
    
    % compute covariance matrix for upper neighbor
    if ~broadband
        fdatHi = filterFGx(data,EEG_in.srate,freqs(i_freq)+neighparams{i_freq}(1,2),neighparams{i_freq}(2,2));
    else
        fdatHi = data;
    end
    %     figure; plot(EEG_in.times,mean(mean(fdatHi(29,:,:),3),1)); title('noise 1')
    fdatHi = reshape( fdatHi(:,time2use,:), EEG_in.nbchan,[] );
    fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
    covHi  = (fdatHi*fdatHi')/sum(time2use);
    % figure; imagesc(covHi); axis square
    
    % combine and regularize noise cov
    ncovm=(covHi+covLo)/2; % combined noise cov matrix
    lambda = regparam * trace(ncovm)/size(ncovm,1);
    % regularisation of noise covariance applied similar to
    % Gulbinaite et al. (2019, NeuroImage)
    % https://doi.org/10.1016/j.neuroimage.2019.116146
    
    % perform generalized eigendecomposition. This is the meat & potatos of RESS
    %     [evecs,evals] = eig(covAt,ncovm); % without regularization
    [evecs,evals] = eig(covAt,ncovm+lambda*eye(size(ncovm))); % with regularization
    [~,comp2plot] = max(diag(evals)); % find maximum component
    evecs = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors (not really necessary, but OK)
    %     figure; plot(sort(diag(evals)),'s-')
    
    % extract components and force sign
    % compmaps = inv(evecs'); % get maps (this is fine for full-rank matrices)
    compmaps = covAt * evecs / (evecs' * covAt * evecs); % this works either way
    [~,idx] = max(abs(compmaps(:,comp2plot))); % find biggest component
    compmaps_n = compmaps * sign(compmaps(idx,comp2plot)); % force to positive sign | how it used to be
    
    % reconstruct RESS component time series
    ress_ts = zeros(EEG_in.pnts,EEG_in.trials);
    for ti=1:EEG_in.trials
        ress_ts(:,ti) = evecs(:,comp2plot)'*squeeze(EEG_in.data(:,:,ti));
    end
    
    % try to subtract component [done by myself] help from Mike Cohen
    noress_ts = zeros(size(EEG_in.data));
    for ti=1:EEG_in.trials
        noress_ts(:,:,ti)  = compmaps(:,setdiff(1:size(evecs,2),comp2plot)) * evecs(:,setdiff(1:size(evecs,2),comp2plot))' * EEG_in.data(:,:,ti);
    end
    
    % save parameters important for exclusion of RESS components later on
    all.compmaps(:,:,i_freq) = compmaps;
    all.comp2plot(i_freq) = comp2plot;
    all.evecs(:,:,i_freq) = evecs;
    
    
    %% calculate SNR-value for amplitude values in FWHM range
    % RESS-component based data
    snr_data(1,:,:)=ress_ts(:,trials2use);
    % filtered data of signal freqs(i_freq) +-peakwidt
    snr_fdatAt = filterFGx(snr_data, EEG_in.srate,freqs(i_freq),peakwidt);
    
    % filtered data of lower and upper neighbor freqs(i_freq) -neighfreq +-peakwidt freqs(i_freq) +neighfreq +-peakwidt
    if ~broadband
        snr_fdatLo = filterFGx(snr_data,EEG_in.srate,freqs(i_freq)+neighparams{i_freq}(1,1),neighparams{i_freq}(2,1));
        snr_fdatHi = filterFGx(snr_data,EEG_in.srate,freqs(i_freq)+neighparams{i_freq}(1,2),neighparams{i_freq}(2,2));
    else
        snr_fdatLo = snr_data;
        snr_fdatHi = snr_data;
    end
    
    % induced and evoked spectra for filtered data
    nfft = ceil( EEG_in.srate/.1 );
    snr_fdatAt_fft_ind = mean(abs( fft(squeeze(snr_fdatAt(1,time2use,:)),nfft,1)*2/sum(time2use) ),2);
    snr_fdatLo_fft_ind = mean(abs( fft(squeeze(snr_fdatLo(1,time2use,:)),nfft,1)*2/sum(time2use) ),2);
    snr_fdatHi_fft_ind = mean(abs( fft(squeeze(snr_fdatHi(1,time2use,:)),nfft,1)*2/sum(time2use) ),2);
    snr_fdatAt_fft_evo = abs( fft(squeeze(mean(snr_fdatAt(1,time2use,:),3)),nfft,2)*2/sum(time2use) );
    snr_fdatLo_fft_evo = abs( fft(squeeze(mean(snr_fdatLo(1,time2use,:),3)),nfft,2)*2/sum(time2use) );
    snr_fdatHi_fft_evo = abs( fft(squeeze(mean(snr_fdatHi(1,time2use,:),3)),nfft,2)*2/sum(time2use) );
    
    % get index for frequencies
    snr_freqs = linspace(0,EEG_in.srate,nfft);
    % calculate some SNR measure based on amplitude of signal divided by amplitude of neighbors
    if ~broadband
        snr_snr_ind = mean(snr_fdatAt_fft_ind(dsearchn(snr_freqs',freqs(i_freq)-peakwidt/2):dsearchn(snr_freqs',freqs(i_freq)+peakwidt/2)))...
            /mean([...
            mean(snr_fdatLo_fft_ind(dsearchn(snr_freqs',freqs(i_freq)+neighparams{i_freq}(1,1)-neighparams{i_freq}(2,1)/2):...
            dsearchn(snr_freqs',freqs(i_freq)+neighparams{i_freq}(1,1)+neighparams{i_freq}(2,1)/2))) ...
            mean(snr_fdatHi_fft_ind(dsearchn(snr_freqs',freqs(i_freq)+neighparams{i_freq}(1,2)-neighparams{i_freq}(2,2)/2):...
            dsearchn(snr_freqs',freqs(i_freq)+neighparams{i_freq}(1,2)+neighparams{i_freq}(2,2)/2)))...
            ]);
        snr_snr_evo = mean(snr_fdatAt_fft_evo(dsearchn(snr_freqs',freqs(i_freq)-peakwidt/2):dsearchn(snr_freqs',freqs(i_freq)+peakwidt/2)))...
            /mean([...
            mean(snr_fdatLo_fft_evo(dsearchn(snr_freqs',freqs(i_freq)+neighparams{i_freq}(1,1)-neighparams{i_freq}(2,1)/2):...
            dsearchn(snr_freqs',freqs(i_freq)+neighparams{i_freq}(1,1)+neighparams{i_freq}(2,1)/2))) ...
            mean(snr_fdatHi_fft_evo(dsearchn(snr_freqs',freqs(i_freq)+neighparams{i_freq}(1,2)-neighparams{i_freq}(2,2)/2):...
            dsearchn(snr_freqs',freqs(i_freq)+neighparams{i_freq}(1,2)+neighparams{i_freq}(2,2)/2)))...
            ]);
    else
        snr_snr_ind = mean(snr_fdatAt_fft_ind(dsearchn(snr_freqs',freqs(i_freq)-peakwidt/2):dsearchn(snr_freqs',freqs(i_freq)+peakwidt/2)))...
            / mean(snr_fdatLo_fft_ind);
        snr_snr_evo = mean(snr_fdatAt_fft_evo(dsearchn(snr_freqs',freqs(i_freq)-peakwidt/2):dsearchn(snr_freqs',freqs(i_freq)+peakwidt/2)))...
            / mean(snr_fdatLo_fft_evo);
    end
    
    
    
    %% plotting for checking/diagnostics
    
    if diagnosticsflag == 1
        figure;
        % topos
        subplot(3,6,[1 2])
        tfix = dsearchn(hz',freqs(i_freq));
        [~,teix]=max(dataXm(:,tfix));
        topoplot(dataX(:,tfix), EEG_in.chanlocs, 'maplimits',[0 max(dataX(:,tfix))],'conv','on', 'numcontour',0,'emarker2',{teix,'o','g',4,1});
        colormap(gca,'hot'); colorbar;
        set(gca,'FontSize',8)
        title(sprintf('induced | %1.3f Hz',freqs(i_freq)),'FontSize',8)
        
        subplot(3,6,[3 4])
        topoplot(dataXm(:,tfix), EEG_in.chanlocs, 'maplimits',[0 max(dataXm(:,tfix))],'conv','on', 'numcontour',0,'emarker2',{teix,'o','g',4,1});
        colormap(gca,'hot'); colorbar;
        set(gca,'FontSize',8)
        title(sprintf('evoked | %1.3f Hz',freqs(i_freq)),'FontSize',8)
        
        subplot(3,6,[5 6])
        map2plot = compmaps_n(:,comp2plot);
        topoplot(map2plot./max(map2plot),EEG_in.chanlocs,'maplimits',[-.7 .7],'numcontour',0,'conv','on','electrodes','off'); colormap(gca,'jet'); colorbar;
        set(gca,'FontSize',8)
        title(sprintf('RESS | %1.3f Hz',freqs(i_freq)),'FontSize',8)
        
        
        % best electrode signal + filtered signal
        subplot(3,6,[7 8])
        [~,teix]=max(dataXm(:,tfix));
        all_trigger=unique(epochvect);
        for i_l = 1:numel(all_trigger)
            plot(EEG_in.times,mean(mean(data(teix,:,epochvect==all_trigger(i_l)),3),1));
            hold on
        end
        ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
        xlim(EEG_in.times([1 end]))
        set(gca,'FontSize',8)
        title(sprintf('signal | %1.3f | best electrode %s', freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
        legend(cellstr(num2str(all_trigger')),'Orientation','horizontal','Location','south','FontSize',6,'box','off')
        
        subplot(3,6,[13 14])
        filtdat = filterFGx(data, EEG_in.srate,freqs(i_freq),peakwidt);
        for i_l = 1:numel(all_trigger)
            plot(EEG_in.times,mean(mean(filtdat(teix,:,epochvect==all_trigger(i_l)),3),1));
            hold on
        end
        ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
        xlim(EEG_in.times([1 end]))
        set(gca,'FontSize',8)
        title(sprintf('filtered signal | %1.3f | best electrode %s', freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
        
        
        % RESS + filtered signal
        subplot(3,6,[9 10])
        all_trigger=unique(epochvect);
        for i_l = 1:numel(all_trigger)
            plot(EEG_in.times,mean(ress_ts(:,epochvect==all_trigger(i_l)),2));
            hold on
        end
        ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
        xlim(EEG_in.times([1 end]))
        set(gca,'FontSize',8)
        title(sprintf('signal | %1.3f | RESS', freqs(i_freq)),'FontSize',8)
        
        
        subplot(3,6,[15 16])
        tdata = reshape(ress_ts,[1 size(ress_ts)]);
        filtdat = filterFGx(tdata, EEG_in.srate,freqs(i_freq),peakwidt);
        for i_l = 1:numel(all_trigger)
            plot(EEG_in.times,mean(mean(filtdat(1,:,epochvect==all_trigger(i_l)),3),1));
            hold on
        end
        ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
        xlim(EEG_in.times([1 end]))
        set(gca,'FontSize',8)
        title(sprintf('filtered signal | %1.3f | RESS', freqs(i_freq)),'FontSize',8)
       
        
        % no RESS + filtered signal
        subplot(3,6,[11 12])
        [~,teix]=max(dataXm(:,tfix));
        all_trigger=unique(epochvect);
        for i_l = 1:numel(all_trigger)
            plot(EEG_in.times,mean(noress_ts(teix,:,epochvect==all_trigger(i_l)),3));
            hold on
        end
        ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
        xlim(EEG_in.times([1 end]))
        set(gca,'FontSize',8)
        title(sprintf('signal | %1.3f | no RESS at %s', freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
        
        
        subplot(3,6,[17 18])
        filtdat = filterFGx(noress_ts, EEG_in.srate,freqs(i_freq),peakwidt);
        for i_l = 1:numel(all_trigger)
            plot(EEG_in.times,mean(mean(filtdat(teix,:,epochvect==all_trigger(i_l)),3),1));
            hold on
        end
        ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
        xlim(EEG_in.times([1 end]))
        set(gca,'FontSize',8)
        title(sprintf('filterd signal | %1.3f | no RESS at %s', freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
        
        
        h.ax1 = axes('Position',[0 .95 1 0.05],'Visible','off');
        if ~broadband
            pl.text = sprintf('\\bf RESS for freq %1.3f | [%1.0f %1.0f]ms | peakwidth %1.2f Hz | neighfreqdist [%1.2f %1.2f] Hz | neighfreqwidth [%1.2f %1.2f] Hz'...
                ,freqs(i_freq),timelims, peakwidt, neighparams{i_freq}(1,:), neighparams{i_freq}(2,:));
        else
            pl.text = sprintf('\\bf RESS for freq %1.3f | [%1.0f %1.0f]ms | peakwidth %1.2f Hz | against broadband'...
                ,freqs(i_freq),timelims, peakwidt);
        end
        text(0.5, 0.5,pl.text,'HorizontalAlignment','center','FontSize',8)
        
        % separate window
        % plotting spectra [best electrode, RESS, noRess]X[evoked induced]
        figure;
        [~,teix]=max(dataXm(:,tfix));
        all_trigger=unique(epochvect);
        pl.xlims = [0 max(freqs)*1.5];
        
        subplot(3,6,[1 2])
        for i_l = 1:numel(all_trigger)
            tdata   = EEG_in.data(teix,time2use,epochvect==all_trigger(i_l));
            tdataX  = mean(abs(fft(tdata,nfft,2)/sum(time2use)).^2,3);
            hz      = linspace(0,EEG_in.srate,nfft);
            plot(hz,tdataX);
            hold on
        end
        xlim(pl.xlims);
        set(gca,'FontSize',8)
        title(sprintf('spectra of induced data | %1.3f Hz | best electrode %s',freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
        legend(cellstr(num2str(all_trigger')),'Orientation','horizontal','Location','north','FontSize',6,'box','off')
        
        
        subplot(3,6,[3 4])
        for i_l = 1:numel(all_trigger)
            tdata   = ress_ts(time2use,epochvect==all_trigger(i_l));
            tdataX  = mean(abs(fft(tdata,nfft,1)/sum(time2use)).^2,2);
            hz      = linspace(0,EEG_in.srate,nfft);
            plot(hz,tdataX);
            hold on
        end
        xlim(pl.xlims);
        set(gca,'FontSize',8)
        title(sprintf('spectra of induced data | %1.3f Hz | RESS',freqs(i_freq)),'FontSize',8)
        
        
        subplot(3,6,[5 6])
        for i_l = 1:numel(all_trigger)
            tdata   = noress_ts(teix,time2use,epochvect==all_trigger(i_l));
            tdataX  = mean(abs(fft(tdata,nfft,2)/sum(time2use)).^2,3);
            hz      = linspace(0,EEG_in.srate,nfft);
            plot(hz,tdataX);
            hold on
        end
        xlim(pl.xlims);
        set(gca,'FontSize',8)
        title(sprintf('spectra of induced data | %1.3f Hz | no RESS best electrode %s',freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
        
        
        subplot(3,6,[7 8])
        for i_l = 1:numel(all_trigger)
            tdata   = EEG_in.data(teix,time2use,epochvect==all_trigger(i_l));
            tdataXm = abs(fft(mean(tdata,3),nfft,2)/sum(time2use)).^2;
            hz      = linspace(0,EEG_in.srate,nfft);
            plot(hz,tdataXm);
            hold on
        end
        xlim(pl.xlims);
        set(gca,'FontSize',8)
        title(sprintf('spectra of evoked data | %1.3f Hz | best electrode %s',freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
        
        
        
        subplot(3,6,[9 10])
        for i_l = 1:numel(all_trigger)
            tdata   = ress_ts(time2use,epochvect==all_trigger(i_l));
            tdataXm = abs(fft(mean(tdata,2),nfft,1)/sum(time2use)).^2;
            hz      = linspace(0,EEG_in.srate,nfft);
            plot(hz,tdataXm);
            hold on
        end
        xlim(pl.xlims);
        set(gca,'FontSize',8)
        title(sprintf('spectra of evoked data | %1.3f Hz | RESS',freqs(i_freq)),'FontSize',8)
        
        
        subplot(3,6,[11 12])
        for i_l = 1:numel(all_trigger)
            tdata   = noress_ts(teix,time2use,epochvect==all_trigger(i_l));
            tdataXm = abs(fft(mean(tdata,3),nfft,2)/sum(time2use)).^2;
            hz      = linspace(0,EEG_in.srate,nfft);
            plot(hz,tdataXm);
            hold on
        end
        xlim(pl.xlims);
        set(gca,'FontSize',8)
        title(sprintf('spectra of evoked data | %1.3f Hz | no RESS best electrode %s',freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
        
        
        subplot(3,6,[13 14 15])
        plot(snr_freqs,snr_fdatAt_fft_ind,'r');
        hold on;
        plot(snr_freqs,snr_fdatLo_fft_ind,'Color',[0 0 0]);
        plot(snr_freqs,snr_fdatHi_fft_ind,'Color',[0.4 0.4 0.4]);
        xlim(pl.xlims);
        set(gca,'FontSize',8)
        title(sprintf('spectra of induced data | RESS | SNR_i_n_d = %1.3f | SNR_e_v_o = %1.3f',snr_snr_ind,snr_snr_evo))
        legend({'induced signal';'induced low-pass noise';'induced high-pass noise'},'FontSize',6,'box','off')
        
        
        subplot(3,6,[16 17 18])
        plot(snr_freqs,snr_fdatAt_fft_evo,'r');
        hold on;
        plot(snr_freqs,snr_fdatLo_fft_evo,'Color',[0 0 0]);
        plot(snr_freqs,snr_fdatHi_fft_evo,'Color',[0.4 0.4 0.4]);
        xlim(pl.xlims);
        set(gca,'FontSize',8)
        title(sprintf('spectra of evoked data | RESS | SNR_i_n_d = %1.3f | SNR_e_v_o = %1.3f',snr_snr_ind,snr_snr_evo))
        legend({'evoked signal';'evoked low-pass noise';'evoked high-pass noise'},'FontSize',6,'box','off')
        
        
        h.ax1 = axes('Position',[0 .95 1 0.05],'Visible','off');
        if ~broadband
            pl.text = sprintf('\\bf RESS for freq %1.3f | [%1.0f %1.0f]ms | peakwidth %1.2f Hz | neighfreqdist [%1.2f %1.2f] Hz | neighfreqwidth [%1.2f %1.2f] Hz'...
                ,freqs(i_freq),timelims, peakwidt, neighparams{i_freq}(1,:), neighparams{i_freq}(2,:));
        else
            pl.text = sprintf('\\bf RESS for freq %1.3f | [%1.0f %1.0f]ms | peakwidth %1.2f Hz | against broadband'...
                ,freqs(i_freq),timelims, peakwidt);
        end
        text(0.5, 0.5,pl.text,'HorizontalAlignment','center','FontSize',8)
    end
    
    %% prepare data for output
    % surrogate EEG lab structure
    RESS{i_freq}.EEG_out = EEG_in;
    RESS{i_freq}.EEG_out = pop_select( RESS{i_freq}.EEG_out,'channel',1);
    RESS{i_freq}.EEG_out = pop_chanedit(RESS{i_freq}.EEG_out, 'changefield',{1 'labels' 'RESS'},...
        'changefield',{1 'theta' ''},'changefield',{1 'radius' ''},'changefield',...
        {1 'X' ''},'changefield',{1 'Y' ''},'changefield',{1 'Z' ''},'changefield',{1 'sph_theta' ''},...
        'changefield',{1 'sph_phi' ''},'changefield',{1 'sph_radius' ''});
    RESS{i_freq}.EEG_out.data(1,:,:) = ress_ts;
    
    % EEG data without RESS component
    RESS{i_freq}.EEG_out_noRESS = EEG_in;
    RESS{i_freq}.EEG_out_noRESS.data = noress_ts;
    
    % EEG out with RESS added as additional channel
    RESS{i_freq}.EEG_out_withRESS = EEG_in;
    RESS{i_freq}.EEG_out_withRESS.data(end+1,:,:) = ress_ts;
    RESS{i_freq}.EEG_out_withRESS.nbchan = EEG_in.nbchan + 1;
    RESS{i_freq}.EEG_out_withRESS.chanlocs(end+1)=RESS{i_freq}.EEG_out_withRESS.chanlocs(end);
    tidx = EEG_in.nbchan + 1;
    tnames = fieldnames(RESS{i_freq}.EEG_out_withRESS.chanlocs);
    tfields = [{''};{'RESS'};repmat({''},numel(tnames)-2,1)];
    for i_field = 1:numel(tnames)
         RESS{i_freq}.EEG_out_withRESS.chanlocs(end).(tnames{i_field})= tfields{i_field};
    end
    
    RESS{i_freq}.EEG_out_withRESS = eeg_checkset(RESS{i_freq}.EEG_out_withRESS);
    
    
    % mean data across conditions
    all_trigger=unique(epochvect);
    for i_l = 1:numel(all_trigger)
        RESS{i_freq}.Mean.data(:,i_l)=mean(ress_ts(:,epochvect==all_trigger(i_l)),2);
        RESS{i_freq}.Mean.condition(i_l)=all_trigger(i_l);
        RESS{i_freq}.Mean.trialnum(i_l)=sum(epochvect==all_trigger(i_l));
    end
    RESS{i_freq}.Mean.times=RESS{i_freq}.EEG_out.times;
    
    % component diagnostics
    RESS{i_freq}.ComponentMaps = compmaps;
    RESS{i_freq}.ComponentMaps_signchanged = compmaps_n;
    RESS{i_freq}.CompNum = comp2plot;
    RESS{i_freq}.eigenvectors = evecs;
    RESS{i_freq}.eigenvalues = diag(evals);
    RESS{i_freq}.Comp_SNR_induced = snr_snr_ind;
    RESS{i_freq}.Comp_SNR_evoked = snr_snr_evo;
    
    % parameters
    RESS{i_freq}.parameters.date = datestr(now);
    RESS{i_freq}.parameters.conds2use = unique(epochvect(trials2use));
    RESS{i_freq}.parameters.trials2use = trials2use;
    RESS{i_freq}.parameters.time = [min(EEG_in.times(time2use)) max(EEG_in.times(time2use))];
    RESS{i_freq}.parameters.peakfreq = freqs(i_freq);
    RESS{i_freq}.parameters.peakwidth = peakwidt;
    RESS{i_freq}.parameters.neighfreqdist = neighfreq;
    RESS{i_freq}.parameters.neighfreqwidth = neighwidt;
    RESS{i_freq}.parameters.neighparams = neighparams;
    RESS{i_freq}.parameters.regparam = regparam;
    RESS{i_freq}.parameters.ress_version =ress_v;
   
end

%% try to remove all RESS components from data...

noress_ts_all = zeros(size(EEG_in.data));
for ti=1:EEG_in.trials
    noress_ts_all(:,:,ti)  = all.compmaps(:,setdiff(1:size(all.evecs(:,:,1),2),all.comp2plot(1)),1)...
        * all.evecs(:,setdiff(1:size(all.evecs(:,:,1),1),all.comp2plot(1)),1)' * EEG_in.data(:,:,ti);
end
for i_freq = 2:numel(freqs)
    for ti=1:EEG_in.trials
        noress_ts_all(:,:,ti)  = all.compmaps(:,setdiff(1:size(all.evecs(:,:,i_freq),2),all.comp2plot(i_freq)),i_freq)...
            * all.evecs(:,setdiff(1:size(all.evecs(:,:,i_freq),2),all.comp2plot(i_freq)),i_freq)' * noress_ts_all(:,:,ti);
    end
end

for i_freq = 1:numel(freqs)
    RESS{i_freq}.EEG_out_noRESS_all = EEG_in;
    RESS{i_freq}.EEG_out_noRESS_all.data = noress_ts_all;
end

% % check plotting
% % no RESS + filtered signal
% figure;
% % best electrode signal + filtered signal
% subplot(3,6,[7 8])
% [~,teix]=max(dataXm(:,tfix));
% all_trigger=unique(epochvect);
% for i_l = 1:numel(all_trigger)
%     plot(EEG_in.times,mean(mean(data(teix,:,epochvect==all_trigger(i_l)),3),1));
%     hold on
% end
% ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
% xlim(EEG_in.times([1 end]))
% set(gca,'FontSize',8)
% title(sprintf('signal | %1.3f | best electrode %s', freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
% legend(cellstr(num2str(all_trigger')),'Orientation','horizontal','Location','south','FontSize',6,'box','off')
% 
% subplot(3,6,[13 14])
% filtdat = filterFGx(data, EEG_in.srate,freqs(i_freq),peakwidt);
% for i_l = 1:numel(all_trigger)
%     plot(EEG_in.times,mean(mean(filtdat(teix,:,epochvect==all_trigger(i_l)),3),1));
%     hold on
% end
% ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
% xlim(EEG_in.times([1 end]))
% set(gca,'FontSize',8)
% title(sprintf('filtered signal | %1.3f | best electrode %s', freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
% 
% 
% % RESS + filtered signal
% subplot(3,6,[9 10])
% all_trigger=unique(epochvect);
% for i_l = 1:numel(all_trigger)
%     plot(EEG_in.times,mean(ress_ts(:,epochvect==all_trigger(i_l)),2));
%     hold on
% end
% ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
% xlim(EEG_in.times([1 end]))
% set(gca,'FontSize',8)
% title(sprintf('signal | %1.3f | RESS', freqs(i_freq)),'FontSize',8)
% 
% 
% subplot(3,6,[15 16])
% tdata = reshape(ress_ts,[1 size(ress_ts)]);
% filtdat = filterFGx(tdata, EEG_in.srate,freqs(i_freq),peakwidt);
% for i_l = 1:numel(all_trigger)
%     plot(EEG_in.times,mean(mean(filtdat(1,:,epochvect==all_trigger(i_l)),3),1));
%     hold on
% end
% ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
% xlim(EEG_in.times([1 end]))
% set(gca,'FontSize',8)
% title(sprintf('filtered signal | %1.3f | RESS', freqs(i_freq)),'FontSize',8)
% 
% % NORESS + filtered signal
% subplot(3,6,[11 12])
% [~,teix]=max(dataXm(:,tfix));
% all_trigger=unique(epochvect);
% for i_l = 1:numel(all_trigger)
%     plot(EEG_in.times,mean(noress_ts_all(teix,:,epochvect==all_trigger(i_l)),3));
%     hold on
%     plot(EEG_in.times,mean(noress_ts(teix,:,epochvect==all_trigger(i_l)),3));
% end
% ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
% xlim(EEG_in.times([1 end]))
% set(gca,'FontSize',8)
% title(sprintf('signal | %1.3f | no RESS at %s', freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
% 
% 
% subplot(3,6,[17 18])
% filtdat = filterFGx(noress_ts_all, EEG_in.srate,freqs(i_freq),peakwidt);
% filtdat2 = filterFGx(noress_ts, EEG_in.srate,freqs(i_freq),peakwidt);
% for i_l = 1:numel(all_trigger)
%     plot(EEG_in.times,mean(mean(filtdat(teix,:,epochvect==all_trigger(i_l)),3),1));
%     hold on
%     plot(EEG_in.times,mean(mean(filtdat2(teix,:,epochvect==all_trigger(i_l)),3),1));
% end
% ylim([-1.1 1.1]*max(cellfun(@(x) max(abs(get(x,'ydata'))),num2cell(get(gca,'Children')))))
% xlim(EEG_in.times([1 end]))
% set(gca,'FontSize',8)
% title(sprintf('filterd signal | %1.3f | no RESS at %s', freqs(i_freq),vararg2str(EEG_in.chanlocs(teix).labels)),'FontSize',8)
% 
% 
% h.ax1 = axes('Position',[0 .95 1 0.05],'Visible','off');
% if ~broadband
%     pl.text = sprintf('\\bf RESS for freq %1.3f | [%1.0f %1.0f]ms | peakwidth %1.2f Hz | neighfreqdist [%1.2f %1.2f] Hz | neighfreqwidth [%1.2f %1.2f] Hz'...
%         ,freqs(i_freq),timelims, peakwidt, neighparams{i_freq}(1,:), neighparams{i_freq}(2,:));
% else
%     pl.text = sprintf('\\bf RESS for freq %1.3f | [%1.0f %1.0f]ms | peakwidth %1.2f Hz | against broadband'...
%         ,freqs(i_freq),timelims, peakwidt);
% end
% text(0.5, 0.5,pl.text,'HorizontalAlignment','center','FontSize',8)

fprintf(1,'...done!\n\n')

end


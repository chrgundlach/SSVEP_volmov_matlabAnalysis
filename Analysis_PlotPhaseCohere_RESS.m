%% plot previously calculated TFA data


clearvars
%% parameters
F.PathIn                = 'E:\work\data\SSVEP_volmov\EEG\PHASECOHERE_RESS_noblinks_CSD';
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 ];
F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks

% Cohere.baseline            = [-2750 -2250];
% Cohere.baseline            = [-2250 -2000];
Cohere.baseline            = [-3500 -3000];
% Cohere.baseline            = [-5000 -4000];
Cohere.baseline            = [-3000 -2750];
Cohere.binnum               = 50;



pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% loop for subjects
for i_sub=1:numel(F.Subjects2Use)
    %%
    % read in TFA data
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_Cohere.mat ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    
    temp.cohere = open(sprintf('%s\\VP%02.0f_exp_Cohere.mat',F.PathIn,F.Subjects2Use(i_sub)));
    
    % preallocate memory
    if i_sub == 1
        Cohere.data_ITC = nan([size(temp.cohere.Cohere.itc) numel(F.Subjects2Use)]);
        Cohere.data_ITC_bc = nan([size(temp.cohere.Cohere.itc) numel(F.Subjects2Use)]);
        Cohere.data_ITCphase = nan([size(temp.cohere.Cohere.itc) numel(F.Subjects2Use)]);
        Cohere.data_ITCphasediff = nan([size(temp.cohere.Cohere.itc)-[1 0] numel(F.Subjects2Use)]);
        Cohere.data_ITCphasediff_bc = nan([size(temp.cohere.Cohere.itc)-[1 0] numel(F.Subjects2Use)]);
        Cohere.data_ERSP_evo = nan([size(temp.cohere.Cohere.ersp_evo) numel(F.Subjects2Use)]);
        Cohere.data_ERSP_ind = nan([size(temp.cohere.Cohere.ersp_ind) numel(F.Subjects2Use)]);
        Cohere.data_SSVEPphaseshift_bc = repmat({[]},1,numel(F.Subjects2Use));
        Cohere.data_SSVEPphaseshift_bc_m = [];
        Cohere.time = temp.cohere.Cohere.times;
        Cohere.frequency = temp.cohere.Cohere.freqs;
        Cohere.electrodes = temp.cohere.Cohere.electrodes;
        Cohere.electrodes_RESS = temp.cohere.Cohere.electrodes_RESS;
        Cohere.RESS_map = nan([length(temp.cohere.Cohere.RESS.map) numel(F.Subjects2Use)]);
        [Cohere.RESS_snr_ind, Cohere.RESS_snr_evo] = deal(nan(1,numel(F.Subjects2Use)));
        Cohere.IndBestElecs = repmat({[]},1,numel(F.Subjects2Use));
        Cohere.data_ITCPhasehist=[];
        Cohere.srate = temp.cohere.Cohere.params.srate;
    end
    
    % subject by subject
    Cohere.data_ERSP_evo(:,:,i_sub) = temp.cohere.Cohere.ersp_evo; % evoked ersp
    Cohere.data_ERSP_evo_bc(:,:,i_sub) = bsxfun(@minus, temp.cohere.Cohere.ersp_evo, ...
        mean(temp.cohere.Cohere.ersp_evo(dsearchn(Cohere.time',Cohere.baseline(1)):dsearchn(Cohere.time',Cohere.baseline(2)),:),1));
    
    Cohere.data_ERSP_ind(:,:,i_sub) = temp.cohere.Cohere.ersp_ind; % induced ersp
    Cohere.data_ERSP_ind_bc(:,:,i_sub) = bsxfun(@minus, temp.cohere.Cohere.ersp_ind, ...
        mean(temp.cohere.Cohere.ersp_ind(dsearchn(Cohere.time',Cohere.baseline(1)):dsearchn(Cohere.time',Cohere.baseline(2)),:),1));
    
    Cohere.data_ITC(:,:,i_sub) = temp.cohere.Cohere.itc; % itc data
    Cohere.data_ITC_bc(:,:,i_sub) = bsxfun(@minus, temp.cohere.Cohere.itc, ...
        mean(temp.cohere.Cohere.itc(dsearchn(Cohere.time',Cohere.baseline(1)):dsearchn(Cohere.time',Cohere.baseline(2)),:),1));
    
    Cohere.data_ITCphase(:,:,i_sub) = ...
         atan(squeeze(mean(cos(temp.cohere.Cohere.itcphase),2))./squeeze(mean(sin(temp.cohere.Cohere.itcphase),2))); % itc phase data
    
    Cohere.data_ITCphasediff(:,:,i_sub) = squeeze(mean(temp.cohere.Cohere.SSVEP_itcphasediff,2)); % itc phase diff data
    Cohere.data_ITCphasediff_bc(:,:,i_sub) = bsxfun(@minus,  Cohere.data_ITCphasediff(:,:,i_sub), ...
        mean(Cohere.data_ITCphasediff(dsearchn(Cohere.time',Cohere.baseline(1)):dsearchn(Cohere.time',Cohere.baseline(2)),:,i_sub),1));
    
    Cohere.IndBestElecs{i_sub}=temp.cohere.Cohere.data_index{2};
    
    Cohere.RESS_map(:,i_sub) = temp.cohere.Cohere.RESS.map(:,end);
    Cohere.RESS_snr_ind(i_sub) = temp.cohere.Cohere.RESS.SNR_ind(end);
    Cohere.RESS_snr_evo(i_sub) = temp.cohere.Cohere.RESS.SNR_evo(end);
    
    % experimental I
    t.tdata=angle(exp(1i*(angle(temp.cohere.Cohere.itc_comp(:,:,1))...
        -repmat(angle(mean(temp.cohere.Cohere.itc_comp(:,:,1),2)),[1,size(temp.cohere.Cohere.itc_comp,2),1]))));
    t.Ndata = nan(numel([-pi:(2*pi/Cohere.binnum):pi])-1,size(t.tdata,1));
    for i_t = 1:size(temp.cohere.Cohere.itcphase,1)
        [t.Ndata(:,i_t),t.edges] = histcounts(t.tdata(i_t,:),[-pi:(2*pi/Cohere.binnum):pi]);
    end
    Cohere.data_ITCPhasehist(:,:,i_sub)=t.Ndata;
    
    % experimental II
%     % extract baseline phase
%     t.ind = dsearchn(Cohere.time',Cohere.baseline(1)):dsearchn(Cohere.time',Cohere.baseline(2));
%     t.data = squeeze(mean(sin(temp.cohere.Cohere.itcphase(t.ind,:,:)),2));
%     t.data = squeeze(sin(temp.cohere.Cohere.itcphase(t.ind,1,:)));
%     % figure; plot(Cohere.time(t.ind),t.data(:,1)); hold on;
%     t.maxval = diff([ min(t.data(:,1)) max(t.data(:,1))])/2;
%     t.sinsig=t.maxval*sin(2*pi*(85/6)*(((0:numel(t.ind)-1))*(1/Cohere.srate)));
%     % plot(Cohere.time(t.ind),t.sinsig);
%     t.convdata = conv(t.sinsig,t.data(:,1)); % figure; plot(t.convdata)
%     t.phaseshift =find(t.convdata==max(t.convdata))-numel(t.sinsig)+1;
%     t.sinsig_sh=t.maxval*sin(2*pi*(85/6)*(((0:numel(t.ind)-1)+t.phaseshift)*(1/Cohere.srate)));
%     % plot(Cohere.time(t.ind),t.sinsig_sh);
%     
%     t.ind2= dsearchn(Cohere.time',Cohere.baseline(1)):numel(Cohere.time);
%     t.data2 = squeeze(mean(sin(temp.cohere.Cohere.itcphase(t.ind2,:,:)),2));
%     t.data2 = squeeze(sin(temp.cohere.Cohere.itcphase(t.ind2,1,:)));
%     t.sinsig_sh2=t.maxval*sin(2*pi*(85/6)*(((0:numel(t.ind2)-1)+t.phaseshift)*(1/Cohere.srate)));
%     % figure; plot(Cohere.time(t.ind2),t.data2(:,1)); hold on;
%     % plot(Cohere.time(t.ind2),t.sinsig_sh2);
%     % figure;
%     % plot(Cohere.time(t.ind2),squeeze(temp.cohere.Cohere.itcphase(t.ind2,1,1))); hold on
%     % plot(Cohere.time(t.ind2),wrapToPi(2*pi*(85/6)*(((0:numel(t.ind2)-1)+t.phaseshift)*(1/Cohere.srate))));
%     figure; plot(Cohere.time(t.ind2),wrapToPi(squeeze(temp.cohere.Cohere.itcphase(t.ind2,1,1))-...
%     wrapToPi(2*pi*(85/6)*(((0:numel(t.ind2)-1)+t.phaseshift)*(1/Cohere.srate)))'));
%     
%     t.fftdata = (fft(t.data,2^12,1))*2/size(t.data,1);
%     t.freqs = ((0:size(t.fftdata,1)-1)/size(t.fftdata,1)) * Cohere.srate;
%     %figure; plot(t.freqs,abs(t.fftdata)); xlim([0 50])
%     %figure; plot(t.freqs,angle(t.fftdata)); xlim([0 50])

    % do loop for trials
    t.ind = dsearchn(Cohere.time',Cohere.baseline(1)):dsearchn(Cohere.time',Cohere.baseline(2));
    t.ind2= dsearchn(Cohere.time',Cohere.baseline(1)):numel(Cohere.time);
    t.resdata = nan([numel(t.ind2) size(temp.cohere.Cohere.itcphase,2),size(temp.cohere.Cohere.itcphase,3)]);
    for i_el = 1:size(temp.cohere.Cohere.itcphase,3)
        for i_tr = 1:size(temp.cohere.Cohere.itcphase,2)
            t.data = squeeze(sin(temp.cohere.Cohere.itcphase(t.ind,i_tr,i_el)));
            t.maxval = diff([ min(t.data) max(t.data)])/2;
            t.sinsig=t.maxval*sin(2*pi*(85/6)*(((0:numel(t.ind)-1))*(1/Cohere.srate)));
            t.convdata = conv(t.sinsig,t.data); % figure; plot(t.convdata)
            t.phaseshift =find(t.convdata==max(t.convdata))-numel(t.data)+1;
            t.sinsig_sh=t.maxval*sin(2*pi*(85/6)*(((0:numel(t.ind)-1)-t.phaseshift)*(1/Cohere.srate)));
            %figure; plot(Cohere.time(t.ind),t.data); hold on; plot(Cohere.time(t.ind),t.sinsig); plot(Cohere.time(t.ind),t.sinsig_sh);
            
            t.phase_dat = squeeze(temp.cohere.Cohere.itcphase(t.ind2,i_tr,i_el));
            t.phase_sh = wrapToPi(2*pi*(85/6)*(((0:numel(t.ind2)-1)-t.phaseshift)*(1/Cohere.srate)))';
            % figure; plot(Cohere.time(t.ind2),t.phase_dat); hold on; plot(Cohere.time(t.ind2),t.phase_sh)
            t.phase_diff = wrapToPi(t.phase_sh-t.phase_dat);
            % figure; plot(Cohere.time(t.ind2), t.phase_diff)
            t.resdata(:,i_tr,i_el)=t.phase_diff;
        end
    end
    Cohere.data_SSVEPphaseshift_bc{1,i_sub}=t.resdata;
    if i_sub ==1; Cohere.data_SSVEPphaseshift_bc_m = nan([size(t.resdata,1) size(t.resdata,3) numel(F.Subjects2Use)]); end;
    Cohere.data_SSVEPphaseshift_bc_m(:,:,i_sub)=squeeze(mean(t.resdata,2));
    Cohere.data_SSVEPphaseshift_bc_m(:,:,i_sub)=squeeze(mean(t.resdata-...
        repmat(mean(t.resdata(1:numel(t.ind),:,:),1),[size(t.resdata,1),1,1]),2));
    % figure; plot(Cohere.time(t.ind2),Cohere.data_SSVEPphaseshift_bc_m)
    % figure; plot(Cohere.time(t.ind2),t.resdata(:,:,1))
    
    
    clear temp
    
end

%% actual plotting RESS data
% plot ERSP values
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
pl.sub2use = 1:numel(F.Subjects2Use);



pl.data = Cohere.data_ERSP_evo(:,1,pl.sub2use);
pl.data_m = mean(pl.data,3);
pl.data_std = std(pl.data,1,3)./sqrt(numel(pl.sub2use)); % standard error

figure;
h.pl=plot(Cohere.time,pl.data_m,'k','LineWidth',1);
hold on;
h.pl_u = plot(Cohere.time,pl.data_m+pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time,pl.data_m-pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('Inter Trial Coherence')
title(sprintf('mean ersp of SSVEPs at RESS component | N = %1.0f',numel(pl.sub2use)))
legend([h.pl h.pl_u],{'mean value';'upper and lower bound | SE'},'Location','SouthOutside','Orientation','horizontal','Box','off')
set(gcf,'Color',[1 1 1])
xlim(pl.xlims)

% plot ERSP values bc
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
% pl.sub2use = F.Subjects2Use;
pl.data = Cohere.data_ERSP_evo_bc(:,1,pl.sub2use);
pl.data_m = mean(pl.data,3);
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(pl.data),0,'dim',2);

figure;
h.pl_vl = plot(Cohere.time,zeros(1,numel(Cohere.time)),'Color',[0.7 0.7 0.7],'LineWidth',0.5);
hold on;
h.pl=plot(Cohere.time,pl.data_m,'k','LineWidth',1);
h.pl_u = plot(Cohere.time,tt.ci(:,2),'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time,tt.ci(:,1),'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('Inter Trial Coherence')
title(sprintf('mean normalized ersp of SSVEPs at RESS component\nN = %1.0f | baseline [%1.0f %1.0f]ms',...
    numel(pl.sub2use),Cohere.baseline))
legend([h.pl h.pl_u],{'mean value';'upper and lower CI | ttest'},'Location','SouthOutside','Orientation','horizontal','Box','off')
xlim(pl.xlims)


% plot ITC values
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
pl.data = Cohere.data_ITC(:,1,pl.sub2use);
pl.data_m = mean(pl.data,3);
pl.data_std = std(pl.data,1,3)./sqrt(numel(pl.sub2use)); % standard error

figure;
h.pl=plot(Cohere.time,pl.data_m,'k','LineWidth',1);
hold on;
h.pl_u = plot(Cohere.time,pl.data_m+pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time,pl.data_m-pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('Inter Trial Coherence')
title(sprintf('mean inter trial coherence of SSVEPs at RESS component | N = %1.0f',numel(pl.sub2use)))
legend([h.pl h.pl_u],{'mean value';'upper and lower bound | SE'},'Location','SouthOutside','Orientation','horizontal','Box','off')
set(gcf,'Color',[1 1 1])
xlim(pl.xlims)

% plot ITC values bc
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
% pl.sub2use = F.Subjects2Use;
pl.data = Cohere.data_ITC_bc(:,1,pl.sub2use);
pl.data_m = mean(pl.data,3);
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(pl.data),0,'dim',2);

figure;
h.pl_vl = plot(Cohere.time,zeros(1,numel(Cohere.time)),'Color',[0.7 0.7 0.7],'LineWidth',0.5);
hold on;
h.pl=plot(Cohere.time,pl.data_m,'k','LineWidth',1);
h.pl_u = plot(Cohere.time,tt.ci(:,2),'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time,tt.ci(:,1),'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('Inter Trial Coherence')
title(sprintf('mean normalized inter trial coherence of SSVEPs at RESS component\nN = %1.0f | baseline [%1.0f %1.0f]ms',...
    numel(pl.sub2use),Cohere.baseline))
legend([h.pl h.pl_u],{'mean value';'upper and lower CI | ttest'},'Location','SouthOutside','Orientation','horizontal','Box','off')
xlim(pl.xlims)

% plot ITC-phase differences
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
% pl.sub2use = F.Subjects2Use;
% pl.data = Cohere.data_ITCphasediff(:,1,pl.sub2use);
pl.data = movmean(Cohere.data_ITCphasediff(:,1,pl.sub2use),10,1);
pl.data_m = mean(pl.data,3);
pl.data_std = std(pl.data,1,3)./sqrt(numel(pl.sub2use)); % standard error

figure;
h.pl=plot(Cohere.time(2:end),pl.data_m,'k','LineWidth',1);
hold on;
h.pl_u = plot(Cohere.time(2:end),pl.data_m+pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time(2:end),pl.data_m-pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('phase lag')
title(sprintf('mean phase lag of SSVEP signals for RESS component | N = %1.0f',numel(pl.sub2use)))
legend([h.pl h.pl_u],{'mean value';'upper and lower bound | SE'},'Location','SouthOutside','Orientation','horizontal','Box','off')
set(gcf,'Color',[1 1 1])
xlim(pl.xlims)

% plot ITC values bc
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
% pl.sub2use = F.Subjects2Use;
% pl.data = Cohere.data_ITCphasediff_bc(:,1,pl.sub2use);
pl.data = movmedian(Cohere.data_ITCphasediff_bc(:,1,pl.sub2use),10,1);
pl.data_m = mean(pl.data,3);
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(pl.data),0,'dim',2);

figure;
h.pl_vl = plot(Cohere.time(2:end),zeros(1,numel(Cohere.time(2:end))),'Color',[0.7 0.7 0.7],'LineWidth',0.5);
hold on;
h.pl=plot(Cohere.time(2:end),pl.data_m,'k','LineWidth',1);
h.pl_u = plot(Cohere.time(2:end),tt.ci(:,2),'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time(2:end),tt.ci(:,1),'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('phase lag')
title(sprintf('mean normalized phase lag of SSVEP signals for RESS component\nN = %1.0f | baseline [%1.0f %1.0f]ms',...
    numel(pl.sub2use),Cohere.baseline))
legend([h.pl h.pl_u],{'mean value';'upper and lower CI | ttest'},'Location','SouthOutside','Orientation','horizontal','Box','off')
xlim(pl.xlims)

%% actual plotting best electrode
% plot ERSP values
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
pl.sub2use = 1:numel(F.Subjects2Use);

pl.data = Cohere.data_ERSP_evo(:,2,pl.sub2use);
pl.data_m = mean(pl.data,3);
pl.data_std = std(pl.data,1,3)./sqrt(numel(pl.sub2use)); % standard error

figure;
h.pl=plot(Cohere.time,pl.data_m,'k','LineWidth',1);
hold on;
h.pl_u = plot(Cohere.time,pl.data_m+pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time,pl.data_m-pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('Inter Trial Coherence')
title(sprintf('mean ersp of SSVEPs at RESS component | N = %1.0f',numel(pl.sub2use)))
legend([h.pl h.pl_u],{'mean value';'upper and lower bound | SE'},'Location','SouthOutside','Orientation','horizontal','Box','off')
set(gcf,'Color',[1 1 1])
xlim(pl.xlims)

% plot ERSP values bc
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
% pl.sub2use = F.Subjects2Use;
pl.data = Cohere.data_ERSP_evo_bc(:,2,pl.sub2use);
pl.data_m = mean(pl.data,3);
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(pl.data),0,'dim',2);

figure;
h.pl_vl = plot(Cohere.time,zeros(1,numel(Cohere.time)),'Color',[0.7 0.7 0.7],'LineWidth',0.5);
hold on;
h.pl=plot(Cohere.time,pl.data_m,'k','LineWidth',1);
h.pl_u = plot(Cohere.time,tt.ci(:,2),'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time,tt.ci(:,1),'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('Inter Trial Coherence')
title(sprintf('mean normalized ersp of SSVEPs at RESS component\nN = %1.0f | baseline [%1.0f %1.0f]ms',...
    numel(pl.sub2use),Cohere.baseline))
legend([h.pl h.pl_u],{'mean value';'upper and lower CI | ttest'},'Location','SouthOutside','Orientation','horizontal','Box','off')
xlim(pl.xlims)


% plot ITC values
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];

pl.data = Cohere.data_ITC(:,2,pl.sub2use);
pl.data_m = mean(pl.data,3);
pl.data_std = std(pl.data,1,3)./sqrt(numel(pl.sub2use)); % standard error

figure;
h.pl=plot(Cohere.time,pl.data_m,'k','LineWidth',1);
hold on;
h.pl_u = plot(Cohere.time,pl.data_m+pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time,pl.data_m-pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('Inter Trial Coherence')
title(sprintf('mean inter trial coherence of SSVEPs at best electrode | N = %1.0f',numel(pl.sub2use)))
legend([h.pl h.pl_u],{'mean value';'upper and lower bound | SE'},'Location','SouthOutside','Orientation','horizontal','Box','off')
set(gcf,'Color',[1 1 1])
xlim(pl.xlims)

% plot ITC values bc
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
% pl.sub2use = F.Subjects2Use;
pl.data = Cohere.data_ITC_bc(:,2,pl.sub2use);
pl.data_m = mean(pl.data,3);
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(pl.data),0,'dim',2);

figure;
h.pl_vl = plot(Cohere.time,zeros(1,numel(Cohere.time)),'Color',[0.7 0.7 0.7],'LineWidth',0.5);
hold on;
h.pl=plot(Cohere.time,pl.data_m,'k','LineWidth',1);
h.pl_u = plot(Cohere.time,tt.ci(:,2),'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time,tt.ci(:,1),'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('Inter Trial Coherence')
title(sprintf('mean normalized inter trial coherence of SSVEPs at best electrode\nN = %1.0f | baseline [%1.0f %1.0f]ms',...
    numel(pl.sub2use),Cohere.baseline))
legend([h.pl h.pl_u],{'mean value';'upper and lower CI | ttest'},'Location','SouthOutside','Orientation','horizontal','Box','off')
xlim(pl.xlims)

% plot ITC-phase differences
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
% pl.sub2use = F.Subjects2Use;
pl.data = Cohere.data_ITCphasediff(:,2,pl.sub2use);
pl.data_m = mean(pl.data,3);
pl.data_std = std(pl.data,1,3)./sqrt(numel(pl.sub2use)); % standard error

figure;
h.pl=plot(Cohere.time(2:end),pl.data_m,'k','LineWidth',1);
hold on;
h.pl_u = plot(Cohere.time(2:end),pl.data_m+pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time(2:end),pl.data_m-pl.data_std,'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('phase lag')
title(sprintf('mean phase lag of SSVEP signals for best electrode | N = %1.0f',numel(pl.sub2use)))
legend([h.pl h.pl_u],{'mean value';'upper and lower bound | SE'},'Location','SouthOutside','Orientation','horizontal','Box','off')
set(gcf,'Color',[1 1 1])
xlim(pl.xlims)

% plot ITC values bc
pl.freq2plot = [85/6 85/6];
pl.xlims = [-5000 5000];
% pl.sub2use = F.Subjects2Use;
pl.data =  Cohere.data_ITCphasediff_bc(:,2,pl.sub2use);
pl.data_m = mean(pl.data,3);
[tt.h,tt.p,tt.ci,tt.stats] = ttest(squeeze(pl.data),0,'dim',2);

figure;
h.pl_vl = plot(Cohere.time(2:end),zeros(1,numel(Cohere.time(2:end))),'Color',[0.7 0.7 0.7],'LineWidth',0.5);
hold on;
h.pl=plot(Cohere.time(2:end),pl.data_m,'k','LineWidth',1);
h.pl_u = plot(Cohere.time(2:end),tt.ci(:,2),'Color',[0.5 0.5 0.5],'LineWidth',1);
h.pl_l = plot(Cohere.time(2:end),tt.ci(:,1),'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel('time in ms')
ylabel('phase lag')
title(sprintf('mean normalized phase lag of SSVEP signals for best electrode\nN = %1.0f | baseline [%1.0f %1.0f]ms',...
    numel(pl.sub2use),Cohere.baseline))
legend([h.pl h.pl_u],{'mean value';'upper and lower CI | ttest'},'Location','SouthOutside','Orientation','horizontal','Box','off')
xlim(pl.xlims)

%% experimental plotting
% first
t.plx=Cohere.time;
t.ply=(-pi:(2*pi/Cohere.binnum):pi)+[(diff(-pi:(2*pi/Cohere.binnum):pi)./2) 0];
t.pldata = sum(Cohere.data_ITCPhasehist,3);
figure; imagesc(t.plx,rad2deg(t.ply(1:end-1)),t.pldata,[0 max(max(t.pldata))]);
set(gca,'ydir','normal'); ylim([-70 70])

% second (diff from average
t.plx=Cohere.time;
t.ply=(-pi:(2*pi/Cohere.binnum):pi)+[(diff(-pi:(2*pi/Cohere.binnum):pi)./2) 0];
t.pldata = bsxfun(@minus, sum(Cohere.data_ITCPhasehist,3), mean(sum(Cohere.data_ITCPhasehist,3),2));
figure; imagesc(t.plx,rad2deg(t.ply(1:end-1)),t.pldata,[-1 1]*max(max(abs(t.pldata))));
set(gca,'ydir','normal'); ylim([-70 70])

%%
pl.freq2plot = [85/6 85/6];
% pl.sub2use = F.Subjects2Use;
pl.data = Cohere.data_ITCphasediff_bc(1,1:end,pl.sub2use);
figure;
h.pl_vl = plot(Cohere.time(2:end),zeros(1,numel(Cohere.time(2:end))),'Color',[0.7 0.7 0.7],'LineWidth',0.5);
hold on;
h.pl=plot(Cohere.time(2:end),squeeze(pl.data),'LineWidth',1);
xlabel('time in ms')
ylabel('Inter Trial Coherence')


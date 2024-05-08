%%%% spectral analysis of pre-TMS time window to get spindle characteristics
%%%% get: spindle frequency, sigma power, sigma amplitude, duration and 1/f level
%%%% and then correlate with trial-by-trial variations in MEP amplitude

%% preliminaries
clc 
clear all; 
close all; 

% fieldtrip
path_ft   = 'C:\Users\siann\Downloads\fieldtrip-20231113\fieldtrip-20231113';
addpath(path_ft);
ft_defaults;

subjects = {'sub-09', 'sub-09', 'sub-06', 'sub-02'}
session = {'ses-exp', 'ses-exp-02', 'ses-exp', 'ses-exp'}
datapath = 'E:\spindle_ppTMS\EEG'
sf = 5000

% define conditions
condition_peak = 'S155'
condition_trough = 'S156'
condition_rising = 'S157'
condition_falling = 'S158'
condition_names = {condition_peak, condition_trough, condition_rising, ...
    condition_falling}

sp_freqs = [12 16]

%% load data
for isub=1:numel(subjects)
    cfg = [];
    cfg.datafile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.eeg']
    cfg.headerfile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.vhdr']
    cfg.continous = 'yes';
    cfg.trialdef.prestim = 2
    cfg.trialdef.poststim = 2 % end segment before TMS pulse to avoid confounding by artifact
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = {condition_peak, condition_trough, condition_falling,...
        condition_rising}
    cfg = ft_definetrial(cfg);
    trial_matrix = cfg.trl


    % read-in data
    cfg.channel = {'C4', 'TP9'};
    cfg.reref = 'yes'
    cfg.refchannel = 'TP9'
    data_raw_c4= ft_preprocessing(cfg);
    
    cfg = [];
    cfg.channel = 'C4'
    data_raw_c4 = ft_selectdata(cfg, data_raw_c4)
    data_raw_c4.label = {'C4'}


    % filter the data in the spindle frequency range
    for itrial=1:length(data_raw_c4.trial)
        data_raw_c4.trial_filtered{itrial} = bandpass(data_raw_c4.trial{itrial}(1,:),sp_freqs,sf)
    end 
    
    % plot all the trials
    % for itrial=1:length(data_raw_c4.trial)
    %     figure;
    %     plot(data_raw_c4.trial{itrial})
    % end 

    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_C4_sp'], "data_raw_c4", '-v7.3')
end

%% TF analysis 
  for isub=1:numel(subjects)
   load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_C4_sp']) % should load datasets with all channels here
    % time-frequency analysis
    % replace data after TMS pulse with zero
    for itrial=1:length(data_raw_c4.trial)
        data_raw_c4.trial{itrial}(1,(2-0.003)*sf:end) = 0
    end 
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning'; % unsure whether hanning or gaussian taper
    cfg.foi          = 1:0.5:35;
    cfg.t_ftimwin    = 5./cfg.foi;  
    cfg.toi          = -1.004:0.02:-0.004; % 1.5 to 1s pre TMS
    cfg.keeptrials   ='yes' % get frequency estimate for every trial
    TFRhann5_all_conds= ft_freqanalysis(cfg, data_raw_c4);  

    cfg              = [];
    cfg.parameter    = 'powspctrm'
    % cfg.baseline     = [-0.5 -0.4] % remove 1/f component from the data 
    cfg.baselinetype = 'absolute';
    cfg.maskstyle    = 'saturation';
    cfg.zlim         = [0 100];
    % cfg.ylim         = [TFRhann5{1}.freq(1,3) TFRhann5{1}.freq(1,end)]
    cfg.channel      = 'C4';
    cfg.interactive  = 'no';
    figure
    ft_singleplotTFR(cfg, TFRhann5_all_conds);
    xlabel('time'); 
    ylabel('frequency');
    title(['time-frequency plot before TMS pulse ' subjects{isub}])

    % check topolplot whether activity in sigma range in really spindle
    % activity 
    % cfg = [];
    % cfg.zlim = [0 50];
    % cfg.xlim = [-1.5 1]; 
    % cfg.ylim = [12 16];
    % % cfg.baseline = [-0.5 -0.4];
    % cfg.baselinetype = 'absolute';
    % cfg.layout = 'acticap-64ch-standard2';
    % figure; ft_topoplotTFR(cfg,TFRhann5_all_conds); colorbar
    % title(['topoplot before TMS pulse'])
        save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_TFRhann5_all_conds'], "TFRhann5_all_conds", '-v7.3')
  end 

%% 
% get RMS value for spindle-frequency filtered signal (sliding window, 200
% ms length)
mean_RMS = []
for isub=1:numel(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_C4_sp']) % load something where data not filtered
    cfg = []
    cfg.latency = [-1.004 -0.004]
    data_raw_c4 = ft_selectdata(cfg, data_raw_c4)
    % filter the data in the spindle frequency range
    for itrial=1:length(data_raw_c4.trial)
        data_raw_c4.trial{itrial} = bandpass(data_raw_c4.trial{itrial}(1,:),sp_freqs,sf)
    end
    Len = 0.2 * sf % window length in samples, 200 ms
    movRMS = dsp.MovingRMS(Len)
    y = cell(1, length(data_raw_c4.trial))
    for itrial=1:numel(data_raw_c4.trial)
        y{itrial} = movRMS(data_raw_c4.trial{itrial}(1,:)) 
        mean_RMS{isub, itrial} = mean(y{itrial})
    end 
 save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, '_mean_RMS_all_subs', subjects{isub}], 'mean_RMS', '-v7.3')
end 

%% get sigma power
 for isub=1:numel(subjects)
   load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_C4_sp'])
   load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub},'_TFRhann5_all_conds'])
   max_sigma = []
   idx_max_sigma = []
   for itrial=1:length(data_raw_c4.trial)
    % find the max power peak in the spindle frequency range 
      [max_sigma(itrial), idx_max_sigma(itrial)] = max([TFRhann5_all_conds.powspctrm(itrial, 1, 23:31, :)], [], 'all') % power for all timepoints, across spindle freq range for C4
   end 
   save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'max_sigma_', subjects{isub}], 'max_sigma', '-v7.3')
   save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'idx_max_sigma_', subjects{isub}], 'idx_max_sigma', '-v7.3')
 end  

%% get spindle frequency
 for isub=1:numel(subjects)
   load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'idx_max_sigma_', subjects{isub}])
   load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub},'_TFRhann5_all_conds'])
   load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_C4_sp'])
    for itrial=1:length(data_raw_c4.trial)
        [I1(itrial),I2(itrial), I3(itrial), I4(itrial)] = ind2sub([1 1 9 51],idx_max_sigma(itrial)) % get idx of max for current trial, C4 channel, originally [1 1 9 51], why?
    end 
    clear idx_max_sigma
    sigma_freqs = TFRhann5_all_conds.freq(23:31)
    max_sigma_freqs = sigma_freqs(I3)
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'max_sigma_freqs_', subjects{isub}], 'max_sigma_freqs', '-v7.3')
 end 

    %% get sigma amplitude
    sigma_amplitude = rms(max_sigma, 1) % rms value for power val associated with peak freq on each trial
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'sigma_amplitude_', subjects{isub}], 'sigma_amplitude', '-v7.3')
    clear max_sigma
    clear max_sigma_freqs
    clear sigma_amplitude
    end


%% get 1/f level
%% how to extract spindle duration??

%% correlate spindle characteristics with trial-by-trial variations in MEP amplitude 
for isub=1:numel(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_MEPs_conds_', subjects{isub}])
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'max_sigma_freqs_', subjects{isub}])
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'max_sigma_', subjects{isub}])
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, '_mean_RMS_all_subs', subjects{isub}])
    % sigma peak freq
    [r_freq{isub}, p_freq{isub}] = corrcoef(max_sigma_freqs, data_MEPs_conds(data_MEPs_conds(:,2)~=159,3)')
    % sigma power
    [r_pow{isub}, p_pow{isub}] = corrcoef(max_sigma, data_MEPs_conds(data_MEPs_conds(:,2)~=159,3)')
    % sigma amplitude / rms
    [r_amp{isub}, p_amp{isub}] = corrcoef(vertcat(mean_RMS{isub,:}), data_MEPs_conds(data_MEPs_conds(:,2)~=159,3))
end 

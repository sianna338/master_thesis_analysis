%%%% post hoc trial sorting for SO activity
%% preliminaires
clc 
clear all; 
close all; 

% fieldtrip
path_ft   = 'C:\Users\siann\Downloads\fieldtrip-20231113\fieldtrip-20231113';
addpath(path_ft);
ft_defaults;

sf = 5000; 

datapath = 'E:\spindle_ppTMS\EEG'

subjects = {'sub-09', 'sub-09', 'sub-06', 'sub-02'}
session = {'ses-exp', 'ses-exp-02', 'ses-exp', 'ses-exp'}
subjects = {'sub-09'}
session = {'ses-exp-02'}

% define conditions
condition_peak = 'S155'
condition_trough = 'S156'
condition_rising = 'S157'
condition_falling = 'S158'

% SO frequency range
so_freqs = [0.16 2]

%% load the data and filter
% look at markers in dataset and segment based on markers
for isub=1:length(subjects)
    cfg = [];
    cfg.datafile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.eeg']
    cfg.headerfile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.vhdr']
    cfg.continous = 'yes';
    cfg.trialdef.prestim = 1.5004
    cfg.trialdef.poststim = -0.004 % end segment before TMS pulse to avoid confounding by artifact
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = {condition_peak, condition_trough, condition_falling,...
        condition_rising}
    cfg = ft_definetrial(cfg);
    trial_matrix = cfg.trl


    % read-in data
    cfg.channel = {'C4', 'TP9'};
    cfg.reref = 'yes'
    cfg.refchannel = 'TP9'
    % cfg.bpfilter = 'yes'
    % cfg.bpfiltord = 2
    % cfg.bpfreq = [10.9 14.9]
    % cfg.bpfiltdir = 'twopass'
    data_raw_c4= ft_preprocessing(cfg);
    
    cfg = [];
    cfg.channel = 'C4'
    data_raw_c4 = ft_selectdata(cfg, data_raw_c4)
    data_raw_c4.label = {'C4'}


    % filter the data in the SO frequency range
    % for itrial=1:length(data_raw_c4.trial)
    %     data_raw_c4.trial{itrial} = bandpass(data_raw_c4.trial{itrial}(1,:),so_freqs,sf)
    % end 
    
    % plot all the trials
    cfg = []
    cfg.demean = 'yes'
    cfg.baseline = [-1.5 -1.004]
    cfg.latency = [-1 -0.004]
    ft_databrowser(cfg, data_raw_c4)
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_C4'], "data_raw_c4", '-v7.3')
end

%% identify slow oscillation events

thresholds = [135.33 141.29 151.15 152.58] % get yasa peak-to_peak amp thresholds

for isub=1:length(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_C4'])
    sign_abs_vector = []
    diff_vector = []
    pos2neg_crossings =  [];
    neg2pos_crossings= [];
    trough_val = []
    peak_val = []
    trough_time = []
    peak_time = []
    for itrial=1:length(data_raw_c4.trial)
        sign_abs_vector{itrial} = sign(data_raw_c4.trial{itrial}(1,:))                                % for each timepoint on each trial, get info whether signal is positive, negative or zero 
        diff_vector{itrial} = diff(sign_abs_vector{itrial})                                           % find differences between abs signal recorded for neighbouring timepoints
        pos2neg_crossings{itrial} =  uint32(find(diff_vector{itrial}<0))
        neg2pos_crossings{itrial} =  uint32(find(diff_vector{itrial}>0))
    % end 
    %   for itrial=1:length(data_raw_c4.trial)
        
     end 
    ptp_value = cell(1, length(data_raw_c4.trial));
    SO_length = cell(1, length(data_raw_c4.trial));
     for itrial=1:length(data_raw_c4.trial)
            for icross=1:length(pos2neg_crossings{itrial}(:))  % for each crossing on the trial
              if ~isempty(pos2neg_crossings{itrial}) && ~isempty(neg2pos_crossings{itrial}) && ...
                length(pos2neg_crossings{itrial}) >= icross && length(neg2pos_crossings{itrial}) >= icross
                if pos2neg_crossings{itrial}(icross) > neg2pos_crossings{itrial}(icross) % make sure pos2neg crossings always preceeds neg2pos % this gives error when poscrossing value empty
                   neg2pos_crossings{itrial}(icross) = [];
                end 
              end 
            end 
     
          if all(pos2neg_crossings{itrial}>1) && ~isempty(neg2pos_crossings{itrial})
            for icross=1:length(pos2neg_crossings{itrial}(:))-1 % for each crossing on the trial
                [trough_val{itrial}(icross), trough_time{itrial}(icross)] = min(data_raw_c4.trial{itrial}(1,pos2neg_crossings{itrial}(icross):neg2pos_crossings{itrial}(icross)))
                [peak_val{itrial}(icross), peak_time{itrial}(icross)] = max(data_raw_c4.trial{itrial}(1,neg2pos_crossings{itrial}(icross):pos2neg_crossings{itrial}(icross+1)))
            end 
          end 
        for icross=1:length(pos2neg_crossings{itrial}(:))-1 % for each crossing on the trial
            trough_time{itrial}(icross) = trough_time{itrial}(icross) + pos2neg_crossings{itrial}(icross) - 1;  % Add offset to get absolute position
            peak_time{itrial}(icross)   = peak_time{itrial}(icross) + neg2pos_crossings{itrial}(icross) - 1;  % Add offset to get absolute position
            ptp_value{itrial}(icross) = abs(trough_val{itrial}(icross))+abs(peak_val{itrial}(icross))      % this seems to be overwirting on each trial, don't know why                              
            SO_length{itrial}(icross)  = double(pos2neg_crossings{itrial}(icross+1) - pos2neg_crossings{itrial}(icross)) / sf; % is this correct?  
        end 
     end  
     save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'ptp_value_', subjects{isub}, '_', session{isub}], "ptp_value", '-v7.3')
     save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'SO_length_', subjects{isub}, '_', session{isub}], "SO_length", '-v7.3')
  end 
        % ptp_value(itrial) = abs(min_val(itrial)) + abs(max_val(itrial));
    %     if (ptp_value(itrial))>thresholds(isub)
    %         trial_threshold_passed = [trial_threshold_passed, itrial]
    %     end 
    %     if ~isempty(trial_threshold_passed)
    %         save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'trial_threshold_passed_', subjects{isub}, '_', session{isub}], "trial_threshold_passed", '-v7.3')
    %     end
    % end 


%% now find for each subject trials where detected SO exceeds amplitude and duration threshold
thresholds = [135.33 141.29 151.15 152.58] % get yasa peak-to_peak amp thresholds
duration_criterion = [0.5 2]
for isub=1:length(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'SO_length_', subjects{isub}, '_', session{isub}])
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'ptp_value_', subjects{isub}, '_', session{isub}])
    trial_threshold_passed = [];
    for itrial=1:length(ptp_value)
        for iSOS=1:length(ptp_value{itrial})
            if any((ptp_value{itrial}(iSOS)>thresholds(isub) & duration_criterion(1)<SO_length{itrial}(iSOS) & ...
                SO_length{itrial}(iSOS)<duration_criterion(2)))
                trial_threshold_passed = [trial_threshold_passed, itrial]
            end 
        end
    end 
    trial_threshold_passed = unique(trial_threshold_passed)
   save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'trial_threshold_passed_', subjects{isub}, '_', session{isub}], "trial_threshold_passed", '-v7.3') 
end 
    
%%
norm = 3 % do we want normalized or raw MEPs?
colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'}; 
overallColor = '#77AC30'
figure; 
hold on; 

all_mean_SO = [];
all_mean_no_SO = [];
data_MEPs_all_no_SO = [];
data_MEPs_all_SO = [];

 for isub=1:length(subjects)
   load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'trial_threshold_passed_', subjects{isub}, '_', session{isub}])
   load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_MEPs_conds_', subjects{isub}])
   data_MEPs_all_SO = [data_MEPs_all_SO; data_MEPs_conds(trial_threshold_passed, norm)]
   data_MEPs_all_no_SO = [data_MEPs_all_no_SO; data_MEPs_conds(setdiff(1:end,trial_threshold_passed),norm)]
   mean_SO{isub} = mean(data_MEPs_conds(trial_threshold_passed,norm))
   mean_no_SO{isub} = mean(data_MEPs_conds(setdiff(1:end,trial_threshold_passed),norm))
   std_SO{isub} = std(data_MEPs_conds(trial_threshold_passed,norm))
   std_no_SO{isub} = std(data_MEPs_conds(setdiff(1:end,trial_threshold_passed),norm))
   counts_threshold_passed{isub} = length(trial_threshold_passed)
   text(0.95, mean_SO{isub}, sprintf('n=%d', counts_threshold_passed{isub}), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top'); 
   handle(isub) = plot([mean_SO{isub}, mean_no_SO{isub}], '-ok', 'MarkerFaceColor',colors{isub},'MarkerSize', 8);
   errorbar([mean_SO{isub}, mean_no_SO{isub}], [std_SO{isub}, std_no_SO{isub}], 'ok', 'MarkerFaceColor', colors{isub}, 'LineWidth', 1, 'Color', colors{isub}); 
 end 

overall_mean = [mean(data_MEPs_all_SO), mean(data_MEPs_all_no_SO)];
overall_std = [std(data_MEPs_all_SO), std(data_MEPs_all_no_SO)];
errorbar([1, 2], overall_mean, overall_std, '-ok', 'MarkerFaceColor', overallColor, ...
    'MarkerSize', 8, 'Color', overallColor, 'LineWidth', 1);
xlim([0.5 2.5]); xticks([1 2]); 
xticklabels({'Slow Oscillation trial', 'no Slow Oscillation Trial'});
ylabel(['MEP amplitude (microvolt)']);
legend({'subject-09 (ses-1)','','subject-09 (ses-2)', '','subject-06', '','subject-02', '', 'Overall Mean'})
grid on; 
hold off; 
meps_nesting_all_subs = struct('subjects', subjects, 'Mean_SO', mean_SO, 'Mean_no_SO', mean_no_SO, 'STD_so', std_SO, 'STD_no_SO', std_no_SO, 'counts_SO_trials', counts_threshold_passed)

%%
[h,p,ci,stats] = ttest2(data_MEPs_all_SO, data_MEPs_all_no_SO)

trials = [4 5 10 13 19 23  38 88 90 99]
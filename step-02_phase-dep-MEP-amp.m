%% preliminaries 
addpath('C:\Users\siann\Downloads\fieldtrip-20231113\fieldtrip-20231113')
ft_defaults

datapath = 'E:\spindle_ppTMS\EEG'
datapath_best = 'E:\spindle_ppTMS\BEST_toolbox'

% define conditions
condition_peak = 'S155'
condition_trough = 'S156'
condition_rising = 'S157'
condition_falling = 'S158'
condition_sp_free = 'S159'

%% 
% concatenate the split recordings from sub-09 
session_sub09 = {'ses-exp-02_01', 'ses-exp-02_02'}; 
for ises = 1:numel(session_sub09)
    filename_sub09{ises} = [datapath, filesep, 'sub-09', filesep, 'ses-exp-02', filesep, 'spindle-ppTMS_', 'sub-09', '_', session_sub09{ises}, '.vhdr']
    hdr{ises} = ft_read_header(filename_sub09{ises}, 'chanindx', 66);
    dat{ises} = ft_read_data(filename_sub09{ises}, 'chanindx', 66);
    evt{ises} = ft_read_event(filename_sub09{ises});
end 

hdr = hdr{1}; 
hdr.nChans = 1
hdr.label = {'FDIr'}
hdr.chantype    = hdr.chantype(1,1)
hdr.chanunit    = hdr.chanunit(1,1)
numSamples1 = size(dat{1}, 2);
dat = cat(2, dat{1}, dat{2});  % concatenate the data along the 2nd dimension

% shift the samples of the second recording to start right after first
% recording
for i=1:length(evt{2})
  evt{2}(i).sample = evt{2}(i).sample + numSamples1;
end

evt = cat(2, evt{1}, evt{2}); % concatenate the events
ft_write_data(([datapath, filesep, 'sub-09', filesep, 'ses-exp-02', filesep, 'spindle-ppTMS_', 'sub-09', '_ses-exp-02.vhdr']), dat, 'header', hdr, 'event', evt); % write data back to disk 
clear hdr
clear dat
clear evt

session_sub06 = {'ses-exp_01', 'ses-exp_02'}; 
for ises = 1:numel(session_sub06)
    filename_sub06{ises} = [datapath, filesep, 'sub-06', filesep, 'ses-exp', filesep, 'spindle-ppTMS_', 'sub-06', '_', session_sub06{ises}, '.vhdr']
    hdr{ises} = ft_read_header(filename_sub06{ises}, 'chanindx', 65);
    dat{ises} = ft_read_data(filename_sub06{ises}, 'chanindx', 65); % only select APB channel because otherwise data too big
    evt{ises} = ft_read_event(filename_sub06{ises});
end 

hdr = hdr{1};   
hdr.nChans = 1
hdr.label = {'APBr'}
hdr.chantype    = hdr.chantype(1,1)
hdr.chanunit    = hdr.chanunit(1,1)
numSamples1 = size(dat{1}, 2);
dat = cat(2, dat{1}, dat{2});  

% shift the samples of the second recording to start right after first
% recording
for i=1:length(evt{2})
  evt{2}(i).sample = evt{2}(i).sample + numSamples1
end

evt = cat(2, evt{1}, evt{2}); % concatenate the events
ft_write_data(([datapath, filesep, 'sub-06', filesep, 'ses-exp', filesep, 'spindle-ppTMS_', 'sub-06', '_ses-exp.vhdr']), dat, 'header', hdr, 'event', evt); % write data back to disk 
clear hdr
clear dat
clear evt

% sub-02
session_sub02 = {'ses-exp_01', 'ses-exp_02', 'ses-exp_03'}; 
for ises = 1:numel(session_sub02)
    filename_sub02{ises} = [datapath, filesep, 'sub-02', filesep, session_sub02{ises}, filesep, 'spindle-ppTMS_', 'sub-02', '_', session_sub02{ises}, '.vhdr']
    hdr{ises} = ft_read_header(filename_sub02{ises}, 'chanindx', 66);
    dat{ises} = ft_read_data(filename_sub02{ises}, 'chanindx', 66); % only select FDI channel because otherwise data too big
    evt{ises} = ft_read_event(filename_sub02{ises});
end 
numSamples1 = size(dat{1}, 2);
numSamples2 = size(dat{2}, 2);
hdr = hdr{1};   
hdr.nChans = 1
hdr.label = {'FDIr'}
hdr.chantype    = hdr.chantype(1,1)
hdr.chanunit    = hdr.chanunit(1,1)
dat = cat(2, dat{1}, dat{2}, dat{3});  

% shift the samples of the second recording to start right after first
% recording
for i=1:length(evt{2})
  evt{2}(i).sample = evt{2}(i).sample + numSamples1
end
for i=1:length(evt{3})
  evt{3}(i).sample = evt{3}(i).sample + numSamples2
end
evt = cat(2, evt{1}, evt{2}, evt{3}); % concatenate the events
ft_write_data(([datapath, filesep, 'sub-02', filesep, 'ses-exp', filesep, 'spindle-ppTMS_', 'sub-02', '_ses-exp.vhdr']), dat, 'header', hdr, 'event', evt); % write data back to disk 
clear hdr
clear dat
clear evt
%%
subjects = {'sub-09', 'sub-09', 'sub-06', 'sub-02'}
session = {'ses-exp', 'ses-exp-02', 'ses-exp', 'ses-exp'}
target_muscle = {'FDIr', 'FDIr', 'APBr', 'FDIr'}
sf = 5000

% look at markers in dataset and segment based on markers
for isub=1:length(subjects)
    cfg = [];
    cfg.datafile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '.eeg']
    cfg.continous = 'yes';
    cfg.channel = target_muscle{isub}
    cfg.trialdef.prestim = 0.05;
    cfg.trialdef.poststim = 0.1;
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = {condition_peak, condition_trough, condition_falling,...
        condition_rising, condition_sp_free}
    cfg = ft_definetrial(cfg);

    % read-in data
    % cfg = [];
    cfg.headerfile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '.vhdr']
    cfg.channel = target_muscle{isub}
    cfg.demean = 'yes';                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                'yes'
    cfg.baselinewindow = [-0.05 -0.005]; 
    data_raw= ft_preprocessing(cfg);
    
    % cfg = []
    % ft_databrowser(cfg, data_raw)
    % extract min and max points for for the time window from 0.015-0.05 after
    % TMS pulse
    cfg=[];
    cfg.latency=[0.015 0.05];
    cfg.channel = target_muscle{isub}
    data_MEPs = ft_selectdata(cfg, data_raw)

    % cfg          = [];
    % cfg.method   = 'trial';
    % cfg.keeptrial = 'nan';
    % data_MEPs_clean = ft_rejectvisual(cfg, data_MEPs);
    % save([datapath, filesep, subjects{isub}, filesep, session{ises}, filesep,'data_', subjects{isub}, '_', session_2{ises}], 'data_MEPs_clean', '-v7.3');

    % extract min and max points for for the time window from 0.015-0.05 after
    % TMS pulse
    for i = 1:numel(data_MEPs.trial)
        min_val = min(data_MEPs.trial{i}(1,:));
        max_val = max(data_MEPs.trial{i}(1,:));
        data_MEPs.mep(i,1) = abs(min_val) + abs(max_val);
    end
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep,'data_MEPs_', subjects{isub}, '_', session{isub}], 'data_MEPs', '-v7.3');
end


%% load the EEG data to check for spindle activity on each stimulation trial
for isub=1:numel(subjects)
    cfg = [];
    cfg.datafile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.eeg']
    cfg.continous = 'yes';
    % cfg.channel = {'all', '-EMG', '-HEOG', '-VEOG', '-ADMr', '-APBr', '-FDIr'}
    cfg.trialdef.prestim = 2;
    cfg.trialdef.poststim = 2;
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = {condition_peak, condition_trough, condition_falling,...
        condition_rising, condition_sp_free}
    cfg = ft_definetrial(cfg);
    trial_matrix = cfg.trl

    cfg.demean = 'yes'
    cfg.baseline = [-0.5 -0.4]; % do baseline correction
    % need to do demeaining to get rid of non-zero DC component in time
    % domain data to not have the TF plot look weird
    %
    % to do: reref to linked mastoids? contralateral mastoid?

    data_tms = ft_preprocessing(cfg);

    % replace data after TMS pulse with zeros
    for itrial=1:length(data_tms.trial)
        data_tms.trial{itrial}(:, ((2-0.004)*sf):end) = 0
    end
    save([datapath, filesep, subjects{isub}, filesep, session{isub},'data_', subjects{isub}, '_', session_2{isub}, 'data_tms'], 'data_tms', '-v7.3')

    % select data for different conditions
    % TFA, frequency dependent window length (5 cycles)
    % multitaper convolution method

    % for itrial = 1:numel(data_c4_pre_tms.trial)
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'wavelet';
    % cfg.width        = 5; % time window depending on frequency, should include 5 cycles
    cfg.taper       = 'hanning';
    cfg.foi          = 1:0.5:20;
    cfg.keeptrials   = 'yes'
    cfg.toi          = -0.5:0.05:-0.01; % should be 0.5 to 0.01 pre TMS
    cfg.pad         = 'nextpow2';
    cfg.polyremoval = 1;
    cfg.width = ceil(cfg.foi * 0.5);     %% number of cycles per frequency for morlet wavelets
    cfg.width(cfg.width < 5) = 5;
    cfg.width(1:8) = [2 3 3 3 3 3 4 4];
    TFRhann5 = ft_freqanalysis(cfg, data_tms);
    % TFRhann5{icond}.freq = round(TFRhann5{icond}.freq*100)/100

    cfg              = [];
    cfg.baselinetype = 'absolute';
    cfg.baseline     = [-0.5 -0.4] % remove 1/f component from the data
    % cfg.maskstyle    = 'saturation';
    % cfg.ylim         = [TFRhann5{1}.freq(1,3) TFRhann5{1}.freq(1,end)]
    cfg.channel      = 'C4';
    cfg.interactive  = 'no';
    figure
    ft_singleplotTFR(cfg, TFRhann5);
    xlabel('time');
    ylabel('frequency');
    title(['time-frequency plot before TMS pulse for channel C4'])

    % plot using imagesc to see whether it looks nicer
    for itrial=1:length(data_tms.trial)
        figure;
        imagesc(TFRhann5.time,TFRhann5.freq,squeeze(TFRhann5.powspctrm(itrial,6,:,:)));axis xy; caxis([-500 30000]);
        xlabel('time');
        ylabel('frequency');
        title(['time-frequency plot before TMS pulse for channel C4, trial:', num2str(itrial), ' condition:', num2str(data_tms.trialinfo(itrial))])
    end

    % topoplot
    % cfg = [];
    % cfg.zlim = [0 10];
    % cfg.xlim = [-1.5 1];
    % cfg.ylim = [1 35];
    % % cfg.baseline = [-0.5 -0.4];
    % cfg.baselinetype = 'absolute';
    % cfg.layout = 'acticap-64ch-standard2';
    % figure; ft_topoplotTFR(cfg,TFRhann5); colorbar
    % title(['topoplot before TMS pulse'])
end


%% this does not work
data_MEPs_conds = [];
MEPs_one_subject = [];
conditions_one_subject = [];
data_normalized_blocks = [];
for isub = 1:numel(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep,'data_MEPs_', subjects{isub},  '_', session{isub}]);
    MEPs_one_subject{isub} = data_MEPs.mep
    conditions_one_subject{isub} = data_MEPs.trialinfo
    data_MEPs_conds = [MEPs_one_subject{isub}, conditions_one_subject{isub}]
    blocks = [1:19:length(data_MEPs.mep)-19] % load data MEPs before and put length here
    for iblock = 1:length(blocks)
        [conditions_counts{iblock},condition_num{iblock}] = groupcounts(data_MEPs_conds(blocks(iblock):(blocks(iblock)+19),2));

        % get mean for each block
        mean_MEP_block{iblock} = nanmean(data_MEPs_conds(blocks(iblock):(blocks(iblock)+19),1))

        % get percent change from mean MEP amplitude for that subject
        data_normalized_blocks{iblock} = ((data_MEPs_conds(blocks(iblock):blocks(iblock)+19,1)-mean_MEP_block{iblock})./mean_MEP_block{iblock})*100 
    end
    % calculate z-scores and remove outliers if needed
    data_MEPs_conds(:,3) = (data_MEPs_conds(:,1)-nanmean(data_MEPs_conds(:,1)))/nanstd(data_MEPs_conds(:,1))
    % threshold = 1.7;
    % non_outliers = (abs(data_MEPs_conds(:,3)) < threshold)
    data_MEPs_conds_normalized{isub} = vertcat(data_normalized_blocks{:})
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, 'data_MEPs_conds_', subjects{isub}], 'data_MEPs_conds', '-v7.3')
end

%% normalize MEPs per block for each subject
data_MEPs_conds_normalized = cell(numel(subjects), 1);

for isub = 1:numel(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_MEPs_', subjects{isub}, '_', session{isub}]);
    MEPs_one_subject = data_MEPs.mep;
    conditions_one_subject = data_MEPs.trialinfo;
    data_MEPs_conds = [MEPs_one_subject, conditions_one_subject]; % store MEP on each trial and condition number
    
    num_blocks = 10;  
    block_size = ceil(length(MEPs_one_subject) / num_blocks); % divide all the trials into blocks 
    data_normalized_blocks = cell(num_blocks, 1); % initialize arrays to store results for each block 
    mean_MEP_block = cell(num_blocks, 1);
    conditions_counts = cell(num_blocks, 1);
    condition_num = cell(num_blocks, 1);
    
    for iblock = 1:num_blocks % for each block
        start_idx = (iblock - 1) * block_size + 1; % ensure each block has right number of trials 
        end_idx = min(iblock * block_size, length(MEPs_one_subject)); % Ensure we don't exceed num of trials
        
        [conditions_counts{iblock}, condition_num{iblock}] = groupcounts(data_MEPs_conds(start_idx:end_idx, 2));
        mean_MEP_block{iblock} = nanmean(data_MEPs_conds(start_idx:end_idx, 1));
        data_normalized_blocks{iblock} = ((data_MEPs_conds(start_idx:end_idx, 1) - mean_MEP_block{iblock}) ./ mean_MEP_block{iblock}) * 100; % get percent change
    end
    
    % Calculate z-scores 
    data_MEPs_conds(:,3) = (MEPs_one_subject - nanmean(MEPs_one_subject)) / nanstd(MEPs_one_subject);
    % non_outliers = abs(z_scores) < threshold; 
    
    data_MEPs_conds_normalized{isub} = vertcat(data_normalized_blocks{:});
    data_MEPs_conds(:,4) = data_MEPs_conds_normalized{isub}
    % Save normalized MEP data
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_MEPs_conds_', subjects{isub}], 'data_MEPs_conds', '-v7.3');
end


%%
% create  overview of MEPs from all subjects
subjects_all = {'sub-09', 'sub-06', 'sub-02'}
data_MEPs = [];
group = []; 
for isub = 1:numel(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep,'data_MEPs_conds_', subjects{isub}])
    data_MEPs{isub} = data_MEPs_conds
    group = [group; repmat(subjects(isub), size(data_MEPs{isub}, 1), 1)];
end 

data_MEPs_all = vertcat(data_MEPs{:}); 
data_MEPs_all(:,5) = categorical(group)
writematrix(data_MEPs_all, [datapath, filesep, 'data_MEPs_4_subs.csv'])

figure; 
boxplot(data_MEPs_all(:,1), group,'Notch','on','Labels',subjects_all)
ylim([0 1500]);
ylabel(['MEP amplitude (micro volt)'])
grid on; 

figure; 
boxplot(data_MEPs_all(:,4), group,'Notch','on','Labels',subjects_all)
% ylim([-100 100]);
ylabel(['MEP amplitude (micro volt) normalized'])
grid on; 
no_trials = groupcounts(group)


%% summary across participants and plot
% Count how many times each condition occurs and then get mean normalized MEP amplitude for each ITI
% [condition_num, ~, idx] = unique(data_MEPs_conds(:,2));

% normalized MEPs
[conditions_counts,condition_num] = groupcounts(data_MEPs_all(:,2)) 
condition_name = ({'peak' 'trough' 'rising' 'falling' 'free'})
condition_num = condition_name
condition_num = categorical(condition_num)
average_MEP = groupsummary(data_MEPs_all(:,4),data_MEPs_all(:,2),"mean")
median_MEP = groupsummary(data_MEPs_all(:,4),data_MEPs_all(:,2),"median")
SD_MEP = groupsummary(data_MEPs_all(:,4),data_MEPs_all(:,2),"std")

MEP_avg_norm = struct('Condition', condition_num, 'Count', conditions_counts, 'Average_MEP', average_MEP, 'SD_MEP', SD_MEP, 'Median_MEP', median_MEP)

% Plot
figure; 
plot(MEP_avg_norm.Condition, MEP_avg_norm.Average_MEP, 'ok', 'MarkerFaceColor','k','MarkerSize', 8); hold on;
errorbar(MEP_avg_norm.Condition, MEP_avg_norm.Average_MEP, MEP_avg_norm.SD_MEP, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]); hold on;
title(['MEP amplitudes normalized (blockwise)'], 'FontWeight', 'bold');
ax = gca
ax.XAxis.Categories={'peak' 'falling' 'trough' 'rising' 'free'};
xlabel(ax,'conditions', 'FontWeight', 'bold');
ylabel(ax,'MEP amplitude (percent chhange from block mean)', 'FontWeight', 'bold');
set(ax, 'FontSize', 12, 'FontName', 'Arial');
grid on; 
set(ax, 'GridLineStyle', '--', 'GridColor', [0.6 0.6 0.6], 'GridAlpha', 0.7);
pbaspect([1.5 1 1])
box on;
colormap(ax, jet)

figure; 
plot(MEP_avg_normalized.Condition, MEP_avg_normalized.Median_MEP, 'ok', 'MarkerFaceColor','k','MarkerSize', 8)
errorbar(MEP_avg_normalized.Condition, MEP_avg_normalized.Median_MEP, SD_MEP, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
title(['Median MEP amplitudes normalized'], 'FontWeight', 'bold');
xlabel('conditions', 'FontWeight', 'bold');
ylabel('median MEP amplitude (microvolt)', 'FontWeight', 'bold');
xticklabels({'peak', '', 'falling', '', 'trough', '', 'rising', '', 'free'})
set(gca, 'FontSize', 12, 'FontName', 'Arial');
grid on; 
set(gca, 'GridLineStyle', '--', 'GridColor', [0.6 0.6 0.6], 'GridAlpha', 0.7);
pbaspect([1.5 1 1]); 
box on;
colormap(jet)

% raw MEPs
[conditions_counts,condition_num] = groupcounts(data_MEPs_all(:,2)) 
condition_name = ({'peak' 'trough' 'rising' 'falling' 'free'})
condition_num = condition_name
condition_num = categorical(condition_num)
average_MEP = groupsummary(data_MEPs_all(:,1),data_MEPs_all(:,2),"mean")
median_MEP = groupsummary(data_MEPs_all(:,1),data_MEPs_all(:,2),"median")
SD_MEP = groupsummary(data_MEPs_all(:,1),data_MEPs_all(:,2),"std")

MEP_avg_raw = struct('Condition', condition_num, 'Count', conditions_counts, 'Average_MEP', average_MEP, 'SD_MEP', SD_MEP, 'Median_MEP', median_MEP)

% Plot
figure; 
plot(MEP_avg_raw.Condition, MEP_avg_raw.Average_MEP, 'ok', 'MarkerFaceColor','k','MarkerSize', 8); hold on;
errorbar(MEP_avg_raw.Condition, MEP_avg_raw.Average_MEP, MEP_avg_raw.SD_MEP, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]); hold on;
title(['MEP amplitudes raw'], 'FontWeight', 'bold');
ax = gca
ax.XAxis.Categories={'peak' 'falling' 'trough' 'rising' 'free'};
xlabel(ax,'conditions', 'FontWeight', 'bold');
ylabel(ax,'MEP amplitude (microvolt)', 'FontWeight', 'bold');
set(ax, 'FontSize', 12, 'FontName', 'Arial');
grid on; 
set(ax, 'GridLineStyle', '--', 'GridColor', [0.6 0.6 0.6], 'GridAlpha', 0.7);
pbaspect([1.5 1 1])
box on;
colormap(ax, jet)


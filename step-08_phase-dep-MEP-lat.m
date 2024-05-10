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

%% load MEP latencies from all subjects and store them with conditions
subjects = {'sub-09', 'sub-09', 'sub-06', 'sub-02'}
session = {'ses-exp', 'ses-exp-02', 'ses-exp', 'ses-exp'}

data_MEPs_conds_lat = cell(numel(subjects), 1);

for isub = 1:numel(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_MEPs_', subjects{isub}, '_', session{isub}]);
    MEPs_lat_one_subject = data_MEPs.latency;
    conditions_one_subject = data_MEPs.trialinfo;
    data_MEPs_conds_lat = [MEPs_lat_one_subject, conditions_one_subject]; % store MEP on each trial and condition number
    
    
    % Calculate z-scores 
    data_MEPs_conds_lat(:,3) = (MEPs_lat_one_subject - nanmean(MEPs_lat_one_subject)) / nanstd(MEPs_lat_one_subject);
    % non_outliers = abs(z_scores) < threshold; 
    
    % Save concatenate MEP data
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_MEPs_conds_lat_', subjects{isub}], 'data_MEPs_conds_lat', '-v7.3');
end

%%
% create  overview of MEPs latencies from all subjects

subjects_all = {'sub-09', 'sub-06', 'sub-02', 'sub-03'}
subjects = {'sub-09', 'sub-09', 'sub-06', 'sub-02', 'sub-03'}
session = {'ses-exp', 'ses-exp-02', 'ses-exp', 'ses-exp', 'ses-exp'}
data_MEPs_lat = [];
group = []; 
for isub = 1:numel(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep,'data_MEPs_conds_lat_', subjects{isub}])
    data_MEPs_lat{isub} = data_MEPs_conds_lat
    group = [group; repmat(subjects(isub), size(data_MEPs_lat{isub}, 1), 1)];
end 

data_MEPs_lat_all = vertcat(data_MEPs_lat{:}); 
data_MEPs_lat_all(:,4) = categorical(group)
writematrix(data_MEPs_lat_all, [datapath, filesep, 'data_MEPs_lat_5subs.csv'])

figure; 
boxplot(data_MEPs_lat_all(:,1), group,'Notch','on','Labels',subjects_all)
ylabel(['MEP latency (s)'])
grid on; 

no_trials = groupcounts(group)


%% summary across participants and plot

% average latency
[conditions_counts,condition_num] = groupcounts(data_MEPs_lat_all(:,2)) 
condition_name = ({'peak' 'trough' 'rising' 'falling' 'free'})
condition_num = condition_name
condition_num = categorical(condition_num)
average_MEP_lat = groupsummary(data_MEPs_lat_all(:,1),data_MEPs_lat_all(:,2),"mean")
median_MEP_lat = groupsummary(data_MEPs_lat_all(:,1),data_MEPs_lat_all(:,2),"median")
SD_MEP_lat = groupsummary(data_MEPs_lat_all(:,1),data_MEPs_lat_all(:,2),"std")

MEP_avg_lat = struct('Condition', condition_num, 'Count', conditions_counts, 'Average_MEP_lat', average_MEP_lat, ...
    'SD_MEP_lat', SD_MEP_lat, 'Median_MEP_lat', median_MEP_lat)

% Plot
figure; 
plot(MEP_avg_lat.Condition, MEP_avg_lat.Average_MEP_lat*1000, 'ok', 'MarkerFaceColor','k','MarkerSize', 8); hold on;
errorbar(MEP_avg_lat.Condition, MEP_avg_lat.Average_MEP_lat*1000, MEP_avg_lat.SD_MEP_lat*1000, 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]); hold on;
title(['MEP latencies'], 'FontWeight', 'bold');
ax = gca
ax.XAxis.Categories={'peak' 'falling' 'trough' 'rising' 'free'};
xlabel(ax,'conditions', 'FontWeight', 'bold');
ylabel(ax,'MEP latency (ms)', 'FontWeight', 'bold');
set(ax, 'FontSize', 12, 'FontName', 'Arial');
grid on; 
set(ax, 'GridLineStyle', '--', 'GridColor', [0.6 0.6 0.6], 'GridAlpha', 0.7);
pbaspect([1.5 1 1])
box on;
colormap(ax, jet)




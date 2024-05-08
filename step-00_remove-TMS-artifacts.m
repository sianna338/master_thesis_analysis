clc 
clear all; 
close all; 

% define conditions
condition_peak = 'S155'
condition_trough = 'S156'
condition_rising = 'S157'
condition_falling = 'S158'
condition_sp_free = 'S159'

%
subjects = {'sub-09', 'sub-09', 'sub-06', 'sub-02'}
session = {'ses-exp', 'ses-exp-02', 'ses-exp', 'ses-exp'}
datapath = 'E:\spindle_ppTMS\EEG'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 1: LOAD IN DATA %%%%%%%%%%%%%%%%%
% look at markers in dataset and segment based on markers
for isub=1:length(subjects)
    cfg = [];
    cfg.datafile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.eeg']
    cfg.headerfile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.vhdr']
    cfg.continous = 'yes';
    cfg.trialdef.prestim = 2.5
    cfg.trialdef.poststim = 2.5
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = {condition_peak, condition_trough, condition_falling,...
        condition_rising, condition_sp_free}
    cfg = ft_definetrial(cfg);
    trial_matrix = cfg.trl
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub} '_trialmatrix'], 'trial_matrix', '-v7.3')


    % read-in data 
    % cfg = []; 
    % cfg.channel = {'C4', 'TP9'};
    % cfg.reref = 'yes'
    % cfg.refchannel = 'TP9'
    data_raw_c4 = ft_preprocessing(cfg);
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_C4'], "data_raw_c4", '-v7.3')
    
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2: Create Event-locked averages %%%%%%%%%%%%%%%%%
    cfg = [];
    cfg.preproc.demean = 'yes';
    cfg.preproc.baselinewindow = [-0.1 -0.001];
    data_tms_avg = ft_timelockanalysis(cfg, data_raw_c4);

    figure;
    plot(data_tms_avg.time, data_tms_avg.avg(2,:)); 
    xlim([-0.1 0.6]);     
    ylim([-40 50]);      
    title(['Channel ' data_tms_avg.label{2}, ' ', subjects{isub}]);
    ylabel('Amplitude (uV)')
    xlabel('Time (s)');
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3: Define artifacts & reject %%%%%%%%%%%%%%%%%
%%
triggers = {condition_peak, condition_trough, condition_falling,...
            condition_rising, condition_sp_free}
for isub=1:length(subjects)
    cfg = [];
    cfg.datafile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.eeg']
    cfg.headerfile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.vhdr']
    cfg.method                  = 'marker'; 
    cfg.prestim                 = .001;     
    cfg.poststim                = .02;    
    cfg.trialdef.eventtype      = 'Stimulus';
    cfg.trialdef.eventvalue     = triggers;
    cfg_ringing = ft_artifact_tms(cfg);  
    
    cfg_artifact = [];
    cfg_artifact.datafile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.eeg']
    cfg_artifact.headerfile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.vhdr']
    cfg_artifact.artfctdef.ringing.artifact = cfg_ringing.artfctdef.tms.artifact; % Add ringing/step response artifact definition

    cfg_artifact.artfctdef.reject = 'partial'; 
    datafile = fullfile(datapath, subjects{isub}, session{isub}, ['data_', subjects{isub}, '_', session{isub}, '_C4.mat']);
    loaded_data = load(datafile, "data_raw_c4");
    cfg_artifact.trl = loaded_data.data_raw_c4.cfg.trl  
    cfg_artifact.artfctdef.minaccepttim = 0.01;
    cfg = ft_rejectartifact(cfg_artifact) % reject parts of trials with TMS artifact, trial number increases
    cfg.channel     = {'all'};
    cfg.reref       = 'yes';
    cfg.refchannel  = {'TP9'};
    cfg.dftfilter = 'yes'
    data_tms_clean  = ft_preprocessing(cfg);
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_ringing_clean'], "data_tms_clean", '-v7.3')
end 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 4: Perform ICA to remove remaining artifacts %%%%%%%%%%%%%%%%%

    cfg = [];
    cfg.demean = 'yes';
    cfg.method = 'fastica';       
    cfg.fastica.approach = 'symm'; 
    cfg.fastica.g = 'gauss';
    
    comp_tms = ft_componentanalysis(cfg, data_tms_clean);
    
    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub},'_comp_tms'], 'comp_tms','-v7.3');
    
    % compute TLA and look at timecourse of each component
    cfg = [];
    comp_tms_avg = ft_timelockanalysis(cfg, comp_tms);
    figure;
    cfg = [];
    cfg.viewmode = 'butterfly';
    ft_databrowser(cfg, comp_tms_avg);
    
    
    % look at topography for components
    comp_tms_2 = comp_tms
    comp_tms_2.unmixing(:,1:64) = comp_tms_2.unmixing(:,1:64)
    comp_tms_2.topolabel = comp_tms_2.topolabel(1:64,:)
    figure;
    cfg           = [];
    cfg.component = [1:69];
    cfg.comment   = 'no';
    cfg.layout    = 'acticap-64ch-standard2.mat'; 
    ft_topoplotIC(cfg, comp_tms_2);
    
    % look at non-timelocked data
    cfg          = [];
    cfg.layout   = 'acticap-64ch-standard2.mat'
    cfg.viewmode = 'component'; % Mode specifically suited to browse through ICA data
    ft_databrowser(cfg, comp_tms_2);
     
    
    % we demeaned data before ICA, so we need to transform raw data agin to
    % component data without demeaning
    cfg          = [];
    cfg.demean   = 'no'; 
    cfg.unmixing = comp_tms.unmixing; 
    cfg.topolabel = comp_tms.topolabel; 
    
    comp_tms          = ft_componentanalysis(cfg, data_tms_clean);
    
    cfg            = [];
    cfg.component  = [ 1 4 39 27];
    cfg.demean     = 'no';
    
    data_components_rejected = ft_rejectcomponent(cfg, comp_tms); % remove components
    
    % look at electrode C4 for cleaned data
    cfg                = [];
    % cfg.vartrllength   = 2;
    cfg.preproc.demean = 'no';
    data_tms_clean_avg = ft_timelockanalysis(cfg, data_component_rejected);

    figure;
    plot(data_tms_clean_avg.time, data_tms_clean_avg.avg(find(contains(data_tms_clean_avg.label,'C4')),:),'b'); % Plot all data
    xlim([-0.1 0.6]); 
    % ylim([-40 50])
    title(['Channel C4']);
    ylabel('Amplitude (uV)')
    xlabel('Time (s)');



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 5: Interpolation %%%%%%%%%%%%%%%%%
for isub=1:length(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_ringing_clean'])
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub} '_trialmatrix'])
    cfg = []
    cfg.trl = trial_matrix
    data_tms_clean = ft_redefinetrial(cfg, data_tms_clean)
    % Interpolate nans using cubic interpolation
    cfg = [];
    cfg.method = 'cubic'; 
    cfg.prewindow = 0.01; % Window prior to segment to use data points for interpolation
    cfg.postwindow = 0.01; % Window after segment to use data points for interpolation
    data_tms_interpolated = ft_interpolatenan(cfg, data_tms_clean); 

    save([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_tms_interpolated'], "data_tms_interpolated", '-v7.3')
    fiff_file  = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_tms_clean_segmented.fif']
    fieldtrip2fiff(fiff_file, data_tms_interpolated)
    % change structure of the data to be one continous recording 
    % num_samples = size(data_tms_clean.trial{1},2)*length(data_tms_clean.trial)
    % trial_matrix_one_trial = [1 num_samples 0 0]
    
    dat = cat(2, data_tms_interpolated.trial{:});  
    datafile = [datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'spindle-ppTMS_', subjects{isub}, '_', session{isub}, '_c4.vhdr']
    hdr = ft_read_header(datafile)
    % clear data_tms_interpolated
    fiff_file  = ([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_tms_clean_unsegmented.vhdr'])
    ft_write_data(fiff_file, dat, 'header', hdr); % write data back to disk 
    % compute the TEP on the cleaned data
    % cfg = [];
    % cfg.preproc.demean = 'yes';
    % cfg.preproc.baselinewindow = [-0.05 -0.001];
    % data_tms_clean_avg = ft_timelockanalysis(cfg, data_tms_interpolated);
    % 
    % figure;
    % plot(data_tms_clean_avg.time, data_tms_clean_avg.avg(find(contains(data_tms_clean_avg.label,'C4')),:),'b'); % Plot all data
    % xlim([-0.1 0.6]); 
    % % ylim([-40 50])
    % title(['Channel C4']);
    % ylabel('Amplitude (uV)')
    % xlabel('Time (s)');
end 


%% exclude spindle-free trials with spindle activity
for isub=1:length(subjects)
    load([datapath, filesep, subjects{isub}, filesep, session{isub}, filesep, 'data_', subjects{isub}, '_', session{isub}, '_tms_interpolated'])
    cfg = []
    cfg.channel = 'C4'
    cfg.demean = 'yes'
    ft_databrowser(cfg, data_tms_interpolated)
end 

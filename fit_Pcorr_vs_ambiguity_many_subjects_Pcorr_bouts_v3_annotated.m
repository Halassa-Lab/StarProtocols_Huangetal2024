%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script analyzes behavioral data from a human discrimination task.
% It identifies periods ("bouts") of high performance based on certain criteria,
% computes performance (P_corr) under different ambiguity levels and conditions,
% and fits psychometric (sigmoid) functions to these data. It then compares
% healthy controls (HCS) and schizophrenia patients (SCZ) on these measures.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up path depending on the system (pc) used.
% The code checks which directory exists to determine the correct file paths 
% for adding necessary functions (plots, utilities) and setting the data path.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if exist('X:\Norman\2CellType\VIP_Arghya\Optotag') ==7 % main pc.
    addpath('Y:\Jonathan')
    addpath('Y:\Jonathan\plots')
    addpath('X:\Norman\2CellType\VIP_Arghya\Optotag')
    drobo_path = 'X:';
elseif exist('Y:\Norman\2CellType\VIP_Arghya\Optotag') ==7 % alt pc/ pc2.
    addpath('Z:\Jonathan')
    addpath('Z:\Jonathan\plots')
    addpath('Y:\Norman\2CellType\VIP_Arghya\Optotag')
    drobo_path = 'Y:';
else % alt pc/ pc3.
    addpath('Y:\Jonathan')
    addpath('Y:\Jonathan\plots')
    addpath('Z:\Norman\2CellType\VIP_Arghya\Optotag')
    drobo_path = 'Z:';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define method to identify performance "bouts":
% Trying several methods - all performance on low conflict trials only.
%%% 1.: From Amanda. Adaptive bout size. Start from trial 1, find trial n as end of bout; start from trial n+1, etc.
%%% 2.: Start from trial of highest performance, expand bout size; then trial of highest performance (outside of previous bouts), rinse and repeat .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method_to_use = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: Set paths, number of pulses, trials, etc.
% n_pulse: total pulses in a trial.
% n_non0_pulse: number of non-zero pulses for evidence manipulation.
% n_trials_task: total number of trials in the task.
% n_trial_else: number of initial trials to ignore (0 here).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_path = strcat(drobo_path,('\Norman\Misc\Human_task\human_amb_data_Anna_Huang\Data\data_20230104\'));

%%% N.B.: Can replace n_pulse & n_non0_pulse to dynamically check Cuefile1-16.
n_pulse = 16;
n_non0_pulse = 9;
n_trials_task = 480;
n_trial_else = 0; %%% Omit no trials.

% Flag to determine whether fitting is done over relative ambiguity or net evidence
flag_fitting_over_rel_amb = 0;

% Flags to control plotting
flag_plot_HCS = 0;
flag_plot_SCZ = 0;
flag_plot_overlayed = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct lists of ambiguity and evidence levels:
% amb_list: ambiguity levels.
% net_evidence_list: net evidence derived from amb_list.
% relative_amb_list: relative ambiguity (amb/ (n_non0 - amb)).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amb_list = [0,2,3,4];
net_evidence_list = n_non0_pulse -2*amb_list;
relative_amb_list = amb_list ./ (n_non0_pulse - amb_list);

% Create smooth lists for plotting fits
net_evidence_smooth_list = 0:0.1:n_non0_pulse;
amb_smooth_list = (n_non0_pulse - net_evidence_smooth_list)/2;
relative_amb_smooth_list = amb_smooth_list ./ (n_non0_pulse - amb_smooth_list);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up parameters for sigmoidal fitting.
% fixed_params_0: Initial fixed parameters for the psychometric function fits.
% For this dataset, constraints are necessary since we have fewer effective points.
%%% N.B.: fitting should not work well as there are basically 3 points (ambiguity level 0 and 2 have P_corr~1), and you need 5 points to fit a sigmoid (and more for psychometric fit function). The current data also only exist for high P_corr so any function will only have a small portion that can be fitted.
%%% With only 3 effective points, need to reduce the number of free parameters and set y_min=0.5, ymax=1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fixed_params_0 = [0.5, NaN, NaN, NaN];
fixed_params_0_congruent = fixed_params_0;
fixed_params_0_congruent(4) = 1;
initial_params_0 = [0.5,0.9,2,0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define anonymous functions for sigmoid fitting:
% sigmoid_fit_function_list: evaluates sigmoid fit over a smooth range.
% sigmoid_fit_function_list_2: evaluates sigmoid fit at discrete net evidence values.
%%% Function fitted, rewrite here for easier plotting.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_fitting_over_rel_amb==0
    sigmoid_fit_function_list = @(Psychometric_params)...
        Psychometric_params(1) + (Psychometric_params(2)-Psychometric_params(1))./(1+10.^((Psychometric_params(3)-net_evidence_smooth_list)*Psychometric_params(4)));
elseif flag_fitting_over_rel_amb==1
    sigmoid_fit_function_list = @(Psychometric_params)...   % Also sigmoid, but xaxis as relative ambiguity,
        Psychometric_params(1) + (Psychometric_params(2)-Psychometric_params(1))./(1+10.^((Psychometric_params(3)-relative_amb_smooth_list)*Psychometric_params(4)));
end

if flag_fitting_over_rel_amb==0
    sigmoid_fit_function_list_2 = @(Psychometric_params)...
        Psychometric_params(1) + (Psychometric_params(2)-Psychometric_params(1))./(1+10.^((Psychometric_params(3)-net_evidence_list)*Psychometric_params(4)));
elseif flag_fitting_over_rel_amb==1
    sigmoid_fit_function_list_2 = @(Psychometric_params)...   % Also sigmoid, but xaxis as relative ambiguity,
        Psychometric_params(1) + (Psychometric_params(2)-Psychometric_params(1))./(1+10.^((Psychometric_params(3)-relative_amb_list)*Psychometric_params(4)));
end

%%% 20210212: Some alternative fit functions... Actually most will not work, need <=3 free params as we only have 4 data points.
sigmoid_fit_function_bias_lapse_list = @(Psychometric_params)...                 % Need to add a parameter 5 = lapse rate = fraction of trials subject choose randomly in sigmfit and fixed_params_0, etc.
    (0.5+Psychometric_params(1)) + (1-Psychometric_params(2) - 0.5-Psychometric_params(1))./(1+10.^((Psychometric_params(3)-net_evidence_smooth_list)*Psychometric_params(4)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify data files for Healthy Control (HCS) and Schizophrenia (SCZ) subjects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initilization

data_list_HCS = dir(strcat(data_path, '10*')); % Healthy Control
data_list_SCZ = dir(strcat(data_path, '20*')); % Schizophrenia patients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate arrays for performance (P_corr) at different ambiguity levels
% and conditions (congruent vs incongruent).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pcorr_amb_list_HCS = NaN(4,length(data_list_HCS));
Pcorr_amb_list_SCZ = NaN(4,length(data_list_SCZ));

Pcorr_amb_congruent_list_HCS = NaN(4,length(data_list_HCS));
Pcorr_amb_congruent_list_SCZ = NaN(4,length(data_list_SCZ));
Pcorr_amb_incongruent_list_HCS = NaN(4,length(data_list_HCS));
Pcorr_amb_incongruent_list_SCZ = NaN(4,length(data_list_SCZ));

param_fit_sigm_list_HCS = NaN(4,length(data_list_HCS));
param_fit_congruent_sigm_list_HCS = NaN(4,length(data_list_HCS));
param_fit_incongruent_sigm_list_HCS = NaN(4,length(data_list_HCS));

param_fit_sigm_list_SCZ = NaN(4,length(data_list_SCZ));
param_fit_congruent_sigm_list_SCZ = NaN(4,length(data_list_SCZ));
param_fit_incongruent_sigm_list_SCZ = NaN(4,length(data_list_SCZ));


%% Loop over data (HCS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over HCS subjects:
% For each subject:
% 1. Read data
% 2. Identify performance bouts using chosen method
% 3. Compute P_corr for different ambiguity levels and conditions.
% 4. Fit a psychometric function to the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-allocate arrays to track, for each subject:
%   (a) the number of incongruent/congruent trials by evidence level (4 levels),
%   (b) the fraction of trials excluded,
%   (c) fraction of "easy" or "hard" trials excluded, etc.
n_incongruent_list_HCS = zeros(4,length(data_list_HCS));
n_congruent_list_HCS = zeros(4,length(data_list_HCS));

P_trials_reject_list_HCS = zeros(1,length(data_list_HCS));
P_easy_trials_reject_list_HCS = zeros(1,length(data_list_HCS));
P_hard_trials_reject_list_HCS = zeros(1,length(data_list_HCS));

% Define thresholds for performance bout finding
Pcorr_smooth_thres = 0.7;
Pcorr_thres = 0.85;
n_smooth_Pcorr = 15;

% Determine the minimal bout size required such that one error does not drop
% performance below Pcorr_thres. The expression:
%    (bout_size_min - 1) / (bout_size_min + 1) > Pcorr_thres
% leads to:
%    bout_size_min > (1 + Pcorr_thres) / (1 - Pcorr_thres).
bout_size_min = ceil((1+Pcorr_thres)/(1-Pcorr_thres)); % Whatever bout size is, should make it such that even with an error on the next trial, Pcorr is still above Pcorr_thres.
%%% In other words, (bout_size_min-1)/(bout_size_min+1)> Pcorr_thres , aka bout_size_min> (1+Pcorr_thres)/(1-Pcorr_thres). 

bout_size_min_2 = 4; % absoulte minimum size of bout. 3?4? .

for i_HCS = 1:length(data_list_HCS)
    % Read subject's data
    data = readtable(strcat(data_path, data_list_HCS(i_HCS).name));

    
    % Extract relevant columns for analysis:
    %   'ambigno' (ambiguous vs. non-ambiguous pulses),
    %   'KeyStimResp_corr' (whether response was correct),
    %   'Congruency' (congruent or incongruent trial type).
    % We skip n_trial_else initial trials if needed.
    n_nonpref_pulse_list = data.ambigno(n_trial_else+1:end);
    response_list = data.KeyStimResp_corr(n_trial_else+1:end);
    congruency_list = data.Congruency(n_trial_else+1:end);
    

    % Optionally, we only consider "low_conflict" trials (e.g., <=2 pulses),
    % setting others to NaN so that we only measure performance in that subset.
    %%%%  if over 2 easiest trials only.    
    response_list_low_conflict = response_list;
    response_list_low_conflict(n_nonpref_pulse_list>2)=NaN;    

    
    %----------------------------------------------------------------------
    % (2) Identify high-performance bouts using chosen method
    %     (method_to_use == 1 or method_to_use == 2)
    %----------------------------------------------------------------------
    if method_to_use==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % METHOD 1: Start from each trial that is correct, then expand or
        % reduce the bout size to remain above Pcorr_thres if possible.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pcorr_bout_start = [];
    Pcorr_bout_end = [];
    i_trial_temp = 1;
    % Move through trials, seeking consecutive blocks of high performance.
    while i_trial_temp < length(response_list_low_conflict)-bout_size_min+1
        % Only attempt to define a bout if the current trial is correct (==1).
        if response_list_low_conflict(i_trial_temp)==1
            i_d_trial_temp = bout_size_min;
            % Compute local performance in the initial prospective bout:
            Pcorr_bout_temp = nanmean(response_list_low_conflict(i_trial_temp:(i_trial_temp+i_d_trial_temp-1)));
            % If performance is above Pcorr_thres, attempt to expand the bout.
            if Pcorr_bout_temp>=Pcorr_thres
                while Pcorr_bout_temp>=Pcorr_thres % If current bout size works, try increase bout size by 1.
                    i_d_trial_temp = i_d_trial_temp + 1;
                    Pcorr_bout_temp = nanmean(response_list_low_conflict(i_trial_temp:(i_trial_temp+i_d_trial_temp-1)));
                    % Stop if we reach the end of available trials.
                    if (i_trial_temp+i_d_trial_temp-1)>=length(response_list_low_conflict)
                        break
                    end
                end
                % We record the final validated block of trials as a high-perf bout.
                Pcorr_bout_start = [Pcorr_bout_start; i_trial_temp];
                Pcorr_bout_end = [Pcorr_bout_end; i_trial_temp+i_d_trial_temp-2]; % Max bout size is i_d_trial_temp-1.
                i_trial_temp = i_trial_temp +i_d_trial_temp-1; %Go to next trial that is outside of bout.
            % If performance is below Pcorr_thres, attempt to shrink the bout.
            elseif Pcorr_bout_temp<Pcorr_thres
                while Pcorr_bout_temp<Pcorr_thres % Try reduce bout size until Pcorr is above threshold, or bout size is too small.
                    i_d_trial_temp = i_d_trial_temp - 1;
                    Pcorr_bout_temp = nanmean(response_list_low_conflict(i_trial_temp:(i_trial_temp+i_d_trial_temp-1)));
                    % If we can re-gain Pcorr_thres by shrinking the block, record it.
                    if Pcorr_bout_temp>=Pcorr_thres
                        Pcorr_bout_start = [Pcorr_bout_start; i_trial_temp];
                        Pcorr_bout_end = [Pcorr_bout_end; i_trial_temp+i_d_trial_temp-2]; % Max bout size is i_d_trial_temp-1.
                        i_trial_temp = i_trial_temp +i_d_trial_temp-1; %Go to next trial that is outside of bout.
                        break
                    end
                    % If the block is too small to maintain Pcorr_thres, give up.
                    if i_d_trial_temp<=bout_size_min_2 % note: must also mean Pcorr_bout_temp<Pcorr_thres.
                        i_trial_temp = i_trial_temp + 1; %Go to next trial that is outside of bout.
                        break
                    end     
                end              
            end
        else
            i_trial_temp = i_trial_temp + 1;
        end
    end
    % Collect all identified bout trials.
    trials_to_use = [];
    for i_bout_temp = 1:length(Pcorr_bout_start)
        trials_to_use = [trials_to_use; [Pcorr_bout_start(i_bout_temp):Pcorr_bout_end(i_bout_temp)]'];
    end
    
    elseif method_to_use==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % METHOD 2: Seek peaks in a smoothed performance vector and expand
        % or contract around that peak to define the "bout".
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pcorr_smooth_list_temp = smoothdata(response_list_low_conflict,'gaussian', n_smooth_Pcorr);    
    Pcorr_bout_max = [];
    Pcorr_bout_start = [];
    Pcorr_bout_end = [];

    % Iterate until no more high-performance bouts can be found
    flag_temp = 1;
    while flag_temp == 1;
        % Identify the maximum smoothed performance in the array.
        [max_temp,i_max_temp] = nanmax(Pcorr_smooth_list_temp);
        % If all are NaN or array is empty, no more bouts can be found.
        if isnan(max_temp)
            flag_temp = 0;
        end
        % Check if the maximum smoothed performance is above threshold.
        if max_temp>=Pcorr_smooth_thres
            Pcorr_bout_max = [Pcorr_bout_max; i_max_temp];
            % For the local max, define an initial block of size "bout_size_min"
            % around i_max_temp, searching for the highest local performance segment.
            Pcorr_list_bout_size_min = NaN(bout_size_min,1);
            for i_bout_temp  = 1:length(Pcorr_list_bout_size_min)
                Pcorr_list_bout_size_min(i_bout_temp) =  nanmean(response_list_low_conflict(nanmax([i_max_temp-bout_size_min+i_bout_temp,1]) : nanmin([i_max_temp-1+i_bout_temp,length(response_list_low_conflict)]) ));
            end
            [~,i_bout_temp] = nanmax(Pcorr_list_bout_size_min);
            i_bout_start_temp = nanmax([i_max_temp-bout_size_min+i_bout_temp,1]);
            i_bout_end_temp = nanmin([i_max_temp-1+i_bout_temp,length(response_list_low_conflict)]);
            i_bout_start_temp_initial = i_bout_start_temp;
            i_bout_end_temp_initial = i_bout_end_temp;
            
            Pcorr_bout_temp = nanmean(response_list_low_conflict(i_bout_start_temp:i_bout_end_temp));
            n_repeat_bout_in_other_bouts = 0;
            % If that local performance is still above Pcorr_thres, we expand
            % the block step by step, ensuring we don't overlap existing bouts.
            if Pcorr_bout_temp>=Pcorr_thres
                while Pcorr_bout_temp >= Pcorr_thres
                    i_bout_start_temp_0 = i_bout_start_temp;
                    i_bout_end_temp_0 = i_bout_end_temp;
    %                 %%%%% If move both start and end of bout per step.
    %                 %%% temp test: if move both the start and end of bout per step, look at and use only the move which leaves the largest Pcorr. Maybe need to do it randomly using MCMC?.
                    %%%%% If move only one of start/end of bout per step, chosen randomly.
                    if rand<=0.5
                        i_bout_start_temp = nanmax(i_bout_start_temp-1,1);
                    else
                        i_bout_end_temp = nanmin(i_bout_end_temp+1,length(response_list_low_conflict));
                    end

                    % Check if this newly enlarged block overlaps any previously stored bout.
                    flag_bout_start_in_other_bouts = 0;
                    for i_bouts = 1:length(Pcorr_bout_start)
                        if (i_bout_start_temp>=Pcorr_bout_start(i_bouts) & i_bout_start_temp<=Pcorr_bout_end(i_bouts)) | (i_bout_end_temp>=Pcorr_bout_start(i_bouts) & i_bout_end_temp<=Pcorr_bout_end(i_bouts))
                            flag_bout_start_in_other_bouts = 1;
                            break
                        end
                    end
                    % If overlap occurs, revert to previous smaller block
                        % and keep track of how many times that happened.
                    if flag_bout_start_in_other_bouts==1
                        i_bout_start_temp = i_bout_start_temp_0;
                        i_bout_end_temp = i_bout_end_temp_0;
                        n_repeat_bout_in_other_bouts = n_repeat_bout_in_other_bouts + 1;
                    end

                    Pcorr_bout_temp = nanmean(response_list_low_conflict(i_bout_start_temp:i_bout_end_temp));
                    % Stop if we've spanned all trials or repeated overlap many times.
                    if (i_bout_end_temp-i_bout_start_temp+1)>=length(response_list_low_conflict)
                        break
                    end
                    if n_repeat_bout_in_other_bouts>=10 % If move results in overlap with other bouts, repeatedly 10 times.
                        break
                    end
                end
                % Store this final expanded high-performance bout.
                Pcorr_bout_start = [Pcorr_bout_start; i_bout_start_temp_0];
                Pcorr_bout_end = [Pcorr_bout_end; i_bout_end_temp_0]; % Max bout size is i_d_trial_temp-1.
                % Zero out that segment in the smoothed array, so it's not re-used.
                Pcorr_smooth_list_temp(i_bout_start_temp_0:i_bout_end_temp_0) = NaN; % Exclude any trials in bout from search of next bout.
               
            % If Pcorr_bout_temp < Pcorr_thres, attempt to shrink the block.    
            elseif Pcorr_bout_temp<Pcorr_thres
                while Pcorr_bout_temp<Pcorr_thres % Try reduce bout size until Pcorr is above threshold, or bout size is too small.
                    if rand<=0.5
                        i_bout_start_temp = nanmax(i_bout_start_temp+1,1);
                    else
                        i_bout_end_temp = nanmin(i_bout_end_temp-1,length(response_list_low_conflict));
                    end
                    %%%%% N.B.: may want to ensure i_max_temp is within the bout.

                    Pcorr_bout_temp = nanmean(response_list_low_conflict(i_bout_start_temp:i_bout_end_temp));
                    if Pcorr_bout_temp>=Pcorr_thres
                        Pcorr_bout_start = [Pcorr_bout_start; i_bout_start_temp];
                        Pcorr_bout_end = [Pcorr_bout_end; i_bout_end_temp]; % Max bout size is i_d_trial_temp-1.
                        Pcorr_smooth_list_temp(i_bout_start_temp_initial:i_bout_end_temp_initial) = NaN; % Exclude any trials in bout from search of next bout.                
                        break
                    end
                    % If the block is now too small, discard.
                    if (i_bout_end_temp-i_bout_start_temp+1)<=bout_size_min_2 % note: must also mean Pcorr_bout_temp<Pcorr_thres.
                        Pcorr_smooth_list_temp(i_bout_start_temp_initial:i_bout_end_temp_initial) = NaN; % Exclude any trials in bout from search of next bout.                
                        break
                    end
                end     
            end
        else
            % If the maximum smoothed performance is below Pcorr_smooth_thres,
            % we are done: no more high-performance bouts to find.
            flag_temp = 0; 
        end    
    end
    % Collect all identified bout trials from method 2.
    trials_to_use = [];
    for i_bout_temp = 1:length(Pcorr_bout_start)
        trials_to_use = [trials_to_use; [Pcorr_bout_start(i_bout_temp):Pcorr_bout_end(i_bout_temp)]'];
    end
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  After finding the trials that belong to high-performance "bouts," we:
%    (1) Restrict the data to those selected trials,
%    (2) Calculate which fraction of the total trials got excluded,
%    (3) Separate trials into congruent vs. incongruent, for further analysis,
%    (4) Compute the mean performance (Pcorr) across different ambiguity levels,
%        as well as the standard errors (SE), and
%    (5) Optionally fit a psychometric function (sigm_fit) for each subject.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Restrict to selected trials (within identified performance bouts)
    trials_to_use = unique(trials_to_use);  
    n_nonpref_pulse_list = n_nonpref_pulse_list(trials_to_use);
    response_list = response_list(trials_to_use);
    congruency_list = congruency_list(trials_to_use);
    
    % Compute proportion of trials excluded
    P_trials_reject_list_HCS(i_HCS) = 1 - length(trials_to_use)/length(response_list_low_conflict );
    P_easy_trials_reject_list_HCS(i_HCS) = 1 - nansum(n_nonpref_pulse_list<=2)/nansum(~isnan(response_list_low_conflict) );
    P_hard_trials_reject_list_HCS(i_HCS) = 1 - nansum(n_nonpref_pulse_list>2)/nansum(isnan(response_list_low_conflict) );
   
    %%% Congruent trials.
    n_nonpref_pulse_list_congruent = n_nonpref_pulse_list(congruency_list=="Congruent");
    response_list_congruent = response_list(congruency_list=="Congruent");
    %%% Incongruent trials.
    n_nonpref_pulse_list_incongruent = n_nonpref_pulse_list(congruency_list=="Incongruent");
    response_list_incongruent = response_list(congruency_list=="Incongruent");
    

    % -----------------------------------------------------------------------
    % For logging the counts of incongruent/congruent trials at each ambiguity level:
    %   (Used later for weighting or verifying balanced design.)
    % -----------------------------------------------------------------------
    n_incongruent_list_HCS(1,i_HCS) = sum(n_nonpref_pulse_list_incongruent==0);
    n_incongruent_list_HCS(2,i_HCS) = sum(n_nonpref_pulse_list_incongruent==2);
    n_incongruent_list_HCS(3,i_HCS) = sum(n_nonpref_pulse_list_incongruent==3);
    n_incongruent_list_HCS(4,i_HCS) = sum(n_nonpref_pulse_list_incongruent==4);
    n_congruent_list_HCS(1,i_HCS) = sum(n_nonpref_pulse_list_congruent==0);
    n_congruent_list_HCS(2,i_HCS) = sum(n_nonpref_pulse_list_congruent==2);
    n_congruent_list_HCS(3,i_HCS) = sum(n_nonpref_pulse_list_congruent==3);
    n_congruent_list_HCS(4,i_HCS) = sum(n_nonpref_pulse_list_congruent==4);


    % -----------------------------------------------------------------------
    % Compute mean performance (Pcorr) and standard error at each ambiguity level
    % for both congruent and incongruent trials.
    %
    %   'amb_list': e.g., [0, 2, 3, 4]
    %   'Pcorr_amb_list_HCS': matrix storing the average performance vs. ambiguity
    % -----------------------------------------------------------------------
    for ind_amb = 1:length(amb_list)
        % Mean performance across all trials (regardless of congruency):
        Pcorr_amb_list_HCS(ind_amb, i_HCS) = nanmean(response_list(n_nonpref_pulse_list==amb_list(ind_amb)));
        % Mean performance for congruent trials only:
        Pcorr_amb_congruent_list_HCS(ind_amb, i_HCS) = nanmean(response_list_congruent(n_nonpref_pulse_list_congruent==amb_list(ind_amb)));
        % Mean performance for incongruent trials only:
        Pcorr_amb_incongruent_list_HCS(ind_amb, i_HCS) = nanmean(response_list_incongruent(n_nonpref_pulse_list_incongruent==amb_list(ind_amb)));
        
        % Corresponding standard errors (SE) of the mean:
        SE_Pcorr_amb_list_HCS(ind_amb, i_HCS) = nanstd(response_list(n_nonpref_pulse_list==amb_list(ind_amb))) / sqrt(size(response_list(n_nonpref_pulse_list==amb_list(ind_amb)),1));
        SE_Pcorr_amb_congruent_list_HCS(ind_amb, i_HCS) = nanstd(response_list_congruent(n_nonpref_pulse_list_congruent==amb_list(ind_amb))) / sqrt(size(response_list_congruent(n_nonpref_pulse_list_congruent==amb_list(ind_amb)),1));
        SE_Pcorr_amb_incongruent_list_HCS(ind_amb, i_HCS) = nanstd(response_list_incongruent(n_nonpref_pulse_list_incongruent==amb_list(ind_amb))) / sqrt(size(response_list_incongruent(n_nonpref_pulse_list_incongruent==amb_list(ind_amb)),1));
    end
    % -----------------------------------------------------------------------
    % If desired, fit a psychometric (sigmoid) function to the performance data
    % at each ambiguity/evidence level. 
    %   'flag_fitting_over_rel_amb' determines whether to fit based on 'net_evidence_list'
    %   or on 'relative_amb_list.'
    %   'sigm_fit' is a custom function or external library function.
    %   'fixed_params_0' / 'fixed_params_0_congruent' are parameter constraints for the fit.
    % -----------------------------------------------------------------------
    if flag_fitting_over_rel_amb==0
        [param_fitted,~,~]              = sigm_fit(net_evidence_list, Pcorr_amb_list_HCS(:,i_HCS)            , fixed_params_0);%, initial_params,0);
        [param_fitted_incongruent,~,~]  = sigm_fit(net_evidence_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS), fixed_params_0);%, initial_params,0);
        [param_fitted_congruent,~,~]    = sigm_fit(net_evidence_list, Pcorr_amb_congruent_list_HCS(:,i_HCS)  , fixed_params_0_congruent);%, initial_params,0);
        
        % Attempt alternative weighting/fitting if parameters are out-of-bounds:   
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit_weight_by_SD(net_evidence_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS), SE_Pcorr_amb_incongruent_list_HCS(:,i_HCS), fixed_params_0);%, initial_params,0);
        end
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit(net_evidence_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS), fixed_params_0);%, initial_params,0);
        end
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit_weight_by_SD(net_evidence_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS), SE_Pcorr_amb_incongruent_list_HCS(:,i_HCS), fixed_params_0, initial_params_0);%, initial_params,0);
        end
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit(net_evidence_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS), fixed_params_0, initial_params_0);%, initial_params,0);
        end

    elseif flag_fitting_over_rel_amb==1
        %%% Fit over relative ambiguity instead
        [param_fitted,~,~]              = sigm_fit(relative_amb_list, Pcorr_amb_list_HCS(:,i_HCS)            , fixed_params_0);%, initial_params,0);
        [param_fitted_congruent,~,~]    = sigm_fit(relative_amb_list, Pcorr_amb_congruent_list_HCS(:,i_HCS)  , fixed_params_0);%, initial_params,0);    
        [param_fitted_incongruent,~,~]  = sigm_fit_weight_by_SD(relative_amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS), SE_Pcorr_amb_incongruent_list_HCS(:,i_HCS), fixed_params_0);%, initial_params,0); 
        
        % Again, attempt alternative weighting/fitting if parameters out-of-range:
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit(relative_amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS), fixed_params_0);%, initial_params,0);
        end    
    end

    % Store the final fitted parameters for this subject (for "all," "congruent," and "incongruent"):
    param_fit_sigm_list_HCS(:,i_HCS) = param_fitted;
    param_fit_congruent_sigm_list_HCS(:,i_HCS) = param_fitted_congruent;
    param_fit_incongruent_sigm_list_HCS(:,i_HCS) = param_fitted_incongruent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Here, we optionally produce diagnostic plots of the psychometric fits
%  for each subject (if flag_plot_HCS==1). We generate a figure with 3x3
%  subplots, each row showing data and fits as a function of:
%    (1) net evidence,
%    (2) absolute ambiguity, and
%    (3) relative ambiguity,
%  and each column grouping "all trials," "incongruent trials," and
%  "congruent trials."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Plots
    if flag_plot_HCS==1
        figure;
        %----------------------------------------------------------------------
        % Subplot (3,3,1): ALL TRIALS vs. NET EVIDENCE
        %   * Plots the fitted curve over 'net_evidence_smooth_list'
        %   * Adds scatter points of actual data vs. 'net_evidence_list'
        %----------------------------------------------------------------------
        subplot(3,3,1)
        hold on;
        plot(net_evidence_smooth_list, sigmoid_fit_function_list(param_fitted));
        scatter(net_evidence_list, Pcorr_amb_list_HCS(:,i_HCS));
        ylabel('Proportion correct');
        xlabel('Net evidence');
        ylim([0.5,1]);
        title('All trials')

        %----------------------------------------------------------------------
        % Subplot (3,3,4): ALL TRIALS vs. ABSOLUTE AMBIGUITY
        %   * Uses 'amb_smooth_list' for the x-axis
        %   * Plots same fitted curve, plus scatter of actual data vs. 'amb_list'
        %----------------------------------------------------------------------
        subplot(3,3,4)
        hold on;
        plot(amb_smooth_list, sigmoid_fit_function_list(param_fitted));
        scatter(amb_list, Pcorr_amb_list_HCS(:,i_HCS));
        ylabel('Proportion correct');
        xlabel('Ambiguity');
        title('All trials')
        ylim([0.5,1]);
    
        %----------------------------------------------------------------------
        % Subplot (3,3,7): ALL TRIALS vs. RELATIVE AMBIGUITY
        %   * Uses 'relative_amb_smooth_list' for the x-axis
        %   * Plots same fitted curve, plus scatter of actual data vs. 'relative_amb_list'
        %----------------------------------------------------------------------
        subplot(3,3,7)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted));
        scatter(relative_amb_list, Pcorr_amb_list_HCS(:,i_HCS));
        ylabel('Proportion correct');
        xlabel('Relative ambiguity');
        title('All trials')
        ylim([0.5,1]);

        %----------------------------------------------------------------------
        % Subplot (3,3,2): INCONGRUENT TRIALS vs. NET EVIDENCE
        %   * We display a separate fit 'param_fitted_incongruent'
        %   * Compare with data from 'Pcorr_amb_incongruent_list_HCS'
        %----------------------------------------------------------------------
        subplot(3,3,2)
        hold on;
        plot(net_evidence_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent));
        scatter(net_evidence_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS));
        ylabel('Proportion correct');
        xlabel('Net evidence');
        title('Incongruent trials')

        %----------------------------------------------------------------------
        % Subplot (3,3,5): INCONGRUENT vs. ABSOLUTE AMBIGUITY
        %----------------------------------------------------------------------
        subplot(3,3,5)
        hold on;
        plot(amb_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent));
        scatter(amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS));
        ylabel('Proportion correct');
        xlabel('Ambiguity');
        title('Incongruent trials')

        %----------------------------------------------------------------------
        % Subplot (3,3,8): INCONGRUENT vs. RELATIVE AMBIGUITY
        %----------------------------------------------------------------------
        subplot(3,3,8)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent));
        scatter(relative_amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS));
        ylabel('Proportion correct');
        xlabel('Relative ambiguity');
        title('Incongruent trials')

        %----------------------------------------------------------------------
        % Subplot (3,3,3): CONGRUENT TRIALS vs. NET EVIDENCE
        %   * Fitted curve is param_fitted_congruent
        %----------------------------------------------------------------------
        subplot(3,3,3)
        hold on;
        plot(net_evidence_smooth_list, sigmoid_fit_function_list(param_fitted_congruent));
        scatter(net_evidence_list, Pcorr_amb_congruent_list_HCS(:,i_HCS));
        ylabel('Proportion correct');
        xlabel('Net evidence');
        title('Congruent trials')
        ylim([0.5,1]);

        %----------------------------------------------------------------------
        % Subplot (3,3,6): CONGRUENT vs. ABSOLUTE AMBIGUITY
        %----------------------------------------------------------------------
        subplot(3,3,6)
        hold on;
        plot(amb_smooth_list, sigmoid_fit_function_list(param_fitted_congruent));
        scatter(amb_list, Pcorr_amb_congruent_list_HCS(:,i_HCS));
        ylabel('Proportion correct');
        xlabel('Ambiguity');
        title('Congruent trials')
        ylim([0.5,1]);

        %----------------------------------------------------------------------
        % Subplot (3,3,9): CONGRUENT vs. RELATIVE AMBIGUITY
        %----------------------------------------------------------------------
        subplot(3,3,9)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted_congruent));
        scatter(relative_amb_list, Pcorr_amb_congruent_list_HCS(:,i_HCS));
        ylabel('Proportion correct');
        xlabel('Relative ambiguity');
        title('Congruent trials')
        ylim([0.5,1]);
    end
end

%% Loop over data (SCZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The same process is repeated for SCZ subjects:
% 1. Identify performance bouts.
% 2. Compute P_corr for different conditions.
% 3. Fit psychometric functions.
% ANNOTATIONS MIRROR THE PREVIOUS SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_incongruent_list_SCZ = zeros(4,length(data_list_SCZ));
n_congruent_list_SCZ = zeros(4,length(data_list_SCZ));

P_trials_reject_list_SCZ = zeros(1,length(data_list_HCS));
P_easy_trials_reject_list_SCZ = zeros(1,length(data_list_HCS));
P_hard_trials_reject_list_SCZ = zeros(1,length(data_list_HCS));

for i_SCZ = 1:length(data_list_SCZ)
    data = readtable(strcat(data_path, data_list_SCZ(i_SCZ).name));

    %%% All trials.
    n_nonpref_pulse_list = data.ambigno(n_trial_else+1:end);
    response_list = data.KeyStimResp_corr(n_trial_else+1:end);
    congruency_list = data.Congruency(n_trial_else+1:end);

    %%%%%% Performance criteria
    %%%%  if over 2 easiest trials only.    
    response_list_low_conflict = response_list;
    response_list_low_conflict(n_nonpref_pulse_list>2)=NaN;    

    % ANNOTATIONS MIRROR THE CORRESPONDING SEGMENTS IN THE HCS ANALYSIS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if method_to_use==1
    Pcorr_bout_start = [];
    Pcorr_bout_end = [];
    i_trial_temp = 1;
    while i_trial_temp < length(response_list_low_conflict)-bout_size_min+1
        if response_list_low_conflict(i_trial_temp)==1
            i_d_trial_temp = bout_size_min;
            Pcorr_bout_temp = nanmean(response_list_low_conflict(i_trial_temp:(i_trial_temp+i_d_trial_temp-1)));
            if Pcorr_bout_temp>=Pcorr_thres
                while Pcorr_bout_temp>=Pcorr_thres % If current bout size works, try increase bout size by 1.
                    i_d_trial_temp = i_d_trial_temp + 1;
                    Pcorr_bout_temp = nanmean(response_list_low_conflict(i_trial_temp:(i_trial_temp+i_d_trial_temp-1)));
                    if (i_trial_temp+i_d_trial_temp-1)>=length(response_list_low_conflict)
                        break
                    end
                end
                Pcorr_bout_start = [Pcorr_bout_start; i_trial_temp];
                Pcorr_bout_end = [Pcorr_bout_end; i_trial_temp+i_d_trial_temp-2]; % Max bout size is i_d_trial_temp-1.
                i_trial_temp = i_trial_temp +i_d_trial_temp-1; %Go to next trial that is outside of bout.
            elseif Pcorr_bout_temp<Pcorr_thres
                while Pcorr_bout_temp<Pcorr_thres % Try reduce bout size until Pcorr is above threshold, or bout size is too small.
%                 while i_d_trial_temp>bout_size_min_2 % Try reduce bout size until Pcorr is above threshold, or bout size is too small.
                    i_d_trial_temp = i_d_trial_temp - 1;
                    Pcorr_bout_temp = nanmean(response_list_low_conflict(i_trial_temp:(i_trial_temp+i_d_trial_temp-1)));
                    if Pcorr_bout_temp>=Pcorr_thres
                        Pcorr_bout_start = [Pcorr_bout_start; i_trial_temp];
                        Pcorr_bout_end = [Pcorr_bout_end; i_trial_temp+i_d_trial_temp-2]; % Max bout size is i_d_trial_temp-1.
                        i_trial_temp = i_trial_temp +i_d_trial_temp-1; %Go to next trial that is outside of bout.
                        break
                    end
                    if i_d_trial_temp<=bout_size_min_2 % note: must also mean Pcorr_bout_temp<Pcorr_thres.
                        i_trial_temp = i_trial_temp + 1; %Go to next trial that is outside of bout.
                        break
                    end
                end    
            end
        else
            i_trial_temp = i_trial_temp + 1;
        end
    end
    
    trials_to_use = [];
    for i_bout_temp = 1:length(Pcorr_bout_start)
        trials_to_use = [trials_to_use; [Pcorr_bout_start(i_bout_temp):Pcorr_bout_end(i_bout_temp)]'];
    end
    
    elseif method_to_use==2
    %%%%% Define performance bout. Method 2.
    Pcorr_smooth_list_temp = smoothdata(response_list_low_conflict,'gaussian', n_smooth_Pcorr);    
    Pcorr_bout_max = [];
    Pcorr_bout_start = [];
    Pcorr_bout_end = [];

    flag_temp = 1;
    while flag_temp == 1;
    
        [max_temp,i_max_temp] = nanmax(Pcorr_smooth_list_temp);
        if isnan(max_temp)
            flag_temp = 0;
        end
        if max_temp>=Pcorr_smooth_thres
            Pcorr_bout_max = [Pcorr_bout_max; i_max_temp];
            %%% First define a bout of size bout_size_min, of highest Pcorr, around i_max_temp.
            Pcorr_list_bout_size_min = NaN(bout_size_min,1);
            for i_bout_temp  = 1:length(Pcorr_list_bout_size_min)
                Pcorr_list_bout_size_min(i_bout_temp) =  nanmean(response_list_low_conflict(nanmax([i_max_temp-bout_size_min+i_bout_temp,1]) : nanmin([i_max_temp-1+i_bout_temp,length(response_list_low_conflict)]) ));
            end
            [~,i_bout_temp] = nanmax(Pcorr_list_bout_size_min);
            i_bout_start_temp = nanmax([i_max_temp-bout_size_min+i_bout_temp,1]);
            i_bout_end_temp = nanmin([i_max_temp-1+i_bout_temp,length(response_list_low_conflict)]);
            i_bout_start_temp_initial = i_bout_start_temp;
            i_bout_end_temp_initial = i_bout_end_temp;
            
            Pcorr_bout_temp = nanmean(response_list_low_conflict(i_bout_start_temp:i_bout_end_temp));
            n_repeat_bout_in_other_bouts = 0;
            
            if Pcorr_bout_temp>=Pcorr_thres
                while Pcorr_bout_temp >= Pcorr_thres
                    i_bout_start_temp_0 = i_bout_start_temp;
                    i_bout_end_temp_0 = i_bout_end_temp;
                    if rand<=0.5
                        i_bout_start_temp = nanmax(i_bout_start_temp-1,1);
                    else
                        i_bout_end_temp = nanmin(i_bout_end_temp+1,length(response_list_low_conflict));
                    end

                    %%%%% N.B.: block if allows bouts to overlap.
                    flag_bout_start_in_other_bouts = 0;
                    for i_bouts = 1:length(Pcorr_bout_start)
                        if (i_bout_start_temp>=Pcorr_bout_start(i_bouts) & i_bout_start_temp<=Pcorr_bout_end(i_bouts)) | (i_bout_end_temp>=Pcorr_bout_start(i_bouts) & i_bout_end_temp<=Pcorr_bout_end(i_bouts))
                            flag_bout_start_in_other_bouts = 1;
                            break
                        end
                    end

                    if flag_bout_start_in_other_bouts==1
                        i_bout_start_temp = i_bout_start_temp_0;
                        i_bout_end_temp = i_bout_end_temp_0;
                        n_repeat_bout_in_other_bouts = n_repeat_bout_in_other_bouts + 1;
                    end

                    Pcorr_bout_temp = nanmean(response_list_low_conflict(i_bout_start_temp:i_bout_end_temp));
                    if (i_bout_end_temp-i_bout_start_temp+1)>=length(response_list_low_conflict)
                        break
                    end
                    if n_repeat_bout_in_other_bouts>=10 % If move results in overlap with other bouts, repeatedly 10 times.
                        break
                    end
                end
                Pcorr_bout_start = [Pcorr_bout_start; i_bout_start_temp_0];
                Pcorr_bout_end = [Pcorr_bout_end; i_bout_end_temp_0]; % Max bout size is i_d_trial_temp-1.
                Pcorr_smooth_list_temp(i_bout_start_temp_0:i_bout_end_temp_0) = NaN; % Exclude any trials in bout from search of next bout.
               
            elseif Pcorr_bout_temp<Pcorr_thres
                while Pcorr_bout_temp<Pcorr_thres % Try reduce bout size until Pcorr is above threshold, or bout size is too small.

                    if rand<=0.5
                        i_bout_start_temp = nanmax(i_bout_start_temp+1,1);
                    else
                        i_bout_end_temp = nanmin(i_bout_end_temp-1,length(response_list_low_conflict));
                    end

                    Pcorr_bout_temp = nanmean(response_list_low_conflict(i_bout_start_temp:i_bout_end_temp));
                    if Pcorr_bout_temp>=Pcorr_thres
                        Pcorr_bout_start = [Pcorr_bout_start; i_bout_start_temp];
                        Pcorr_bout_end = [Pcorr_bout_end; i_bout_end_temp]; % Max bout size is i_d_trial_temp-1.
                        Pcorr_smooth_list_temp(i_bout_start_temp_initial:i_bout_end_temp_initial) = NaN; % Exclude any trials in bout from search of next bout.                
                        break
                    end
                    if (i_bout_end_temp-i_bout_start_temp+1)<=bout_size_min_2 % note: must also mean Pcorr_bout_temp<Pcorr_thres.
                        Pcorr_smooth_list_temp(i_bout_start_temp_initial:i_bout_end_temp_initial) = NaN; % Exclude any trials in bout from search of next bout.                
                        break
                    end      
                end  
            end
        else
            flag_temp = 0; % Stop defining new bouts.
        end    
    end
    
    trials_to_use = [];
    for i_bout_temp = 1:length(Pcorr_bout_start)
        trials_to_use = [trials_to_use; [Pcorr_bout_start(i_bout_temp):Pcorr_bout_end(i_bout_temp)]'];
    end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    trials_to_use = unique(trials_to_use);
    n_nonpref_pulse_list = n_nonpref_pulse_list(trials_to_use);
    response_list = response_list(trials_to_use);
    congruency_list = congruency_list(trials_to_use);
    
    P_trials_reject_list_SCZ(i_SCZ) = 1 - length(trials_to_use)/length(response_list_low_conflict );
    P_easy_trials_reject_list_SCZ(i_SCZ) = 1 - nansum(n_nonpref_pulse_list<=2)/nansum(~isnan(response_list_low_conflict) );
    P_hard_trials_reject_list_SCZ(i_SCZ) = 1 - nansum(n_nonpref_pulse_list>2)/nansum(isnan(response_list_low_conflict) );
    
    %%% Congruent trials.
    n_nonpref_pulse_list_congruent = n_nonpref_pulse_list(congruency_list=="Congruent");
    response_list_congruent = response_list(congruency_list=="Congruent");
    %%% Incongruent trials.
    n_nonpref_pulse_list_incongruent = n_nonpref_pulse_list(congruency_list=="Incongruent");
    response_list_incongruent = response_list(congruency_list=="Incongruent");

    n_incongruent_list_SCZ(1,i_SCZ) = sum(n_nonpref_pulse_list_incongruent==0);
    n_incongruent_list_SCZ(2,i_SCZ) = sum(n_nonpref_pulse_list_incongruent==2);
    n_incongruent_list_SCZ(3,i_SCZ) = sum(n_nonpref_pulse_list_incongruent==3);
    n_incongruent_list_SCZ(4,i_SCZ) = sum(n_nonpref_pulse_list_incongruent==4);
    n_congruent_list_SCZ(1,i_SCZ) = sum(n_nonpref_pulse_list_congruent==0);
    n_congruent_list_SCZ(2,i_SCZ) = sum(n_nonpref_pulse_list_congruent==2);
    n_congruent_list_SCZ(3,i_SCZ) = sum(n_nonpref_pulse_list_congruent==3);
    n_congruent_list_SCZ(4,i_SCZ) = sum(n_nonpref_pulse_list_congruent==4);
    
    for ind_amb = 1:length(amb_list)
        Pcorr_amb_list_SCZ(ind_amb, i_SCZ) = nanmean(response_list(n_nonpref_pulse_list==amb_list(ind_amb)));
        Pcorr_amb_congruent_list_SCZ(ind_amb, i_SCZ) = nanmean(response_list_congruent(n_nonpref_pulse_list_congruent==amb_list(ind_amb)));
        Pcorr_amb_incongruent_list_SCZ(ind_amb, i_SCZ) = nanmean(response_list_incongruent(n_nonpref_pulse_list_incongruent==amb_list(ind_amb)));
                
        SE_Pcorr_amb_list_SCZ(ind_amb, i_SCZ) = nanstd(response_list(n_nonpref_pulse_list==amb_list(ind_amb))) / sqrt(size(response_list(n_nonpref_pulse_list==amb_list(ind_amb)),1));
        SE_Pcorr_amb_congruent_list_SCZ(ind_amb, i_SCZ) = nanstd(response_list_congruent(n_nonpref_pulse_list_congruent==amb_list(ind_amb))) / sqrt(size(response_list_congruent(n_nonpref_pulse_list_congruent==amb_list(ind_amb)),1));
        SE_Pcorr_amb_incongruent_list_SCZ(ind_amb, i_SCZ) = nanstd(response_list_incongruent(n_nonpref_pulse_list_incongruent==amb_list(ind_amb))) / sqrt(size(response_list_incongruent(n_nonpref_pulse_list_incongruent==amb_list(ind_amb)),1));
    end

    if flag_fitting_over_rel_amb==0
        [param_fitted,~,~]              = sigm_fit(net_evidence_list, Pcorr_amb_list_SCZ(:,i_SCZ)            , fixed_params_0);%, initial_params,0);
        [param_fitted_incongruent,~,~]  = sigm_fit(net_evidence_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), fixed_params_0);%, initial_params,0);
        [param_fitted_congruent,~,~]    = sigm_fit(net_evidence_list, Pcorr_amb_congruent_list_SCZ(:,i_SCZ)  , fixed_params_0_congruent);%, initial_params,0);
     
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit_weight_by_SD(net_evidence_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), SE_Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), fixed_params_0);%, initial_params,0);
        end
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit(net_evidence_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), fixed_params_0);%, initial_params,0);
        end
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit_weight_by_SD(net_evidence_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), SE_Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), fixed_params_0, initial_params_0);%, initial_params,0);
        end
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit(net_evidence_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), fixed_params_0, initial_params_0);%, initial_params,0);
        end
        
    elseif flag_fitting_over_rel_amb==1
        %%% Fit over relative ambiguity instead
        [param_fitted,~,~]              = sigm_fit(relative_amb_list, Pcorr_amb_list_SCZ(:,i_SCZ)            , fixed_params_0);%, initial_params,0);
        [param_fitted_congruent,~,~]    = sigm_fit(relative_amb_list, Pcorr_amb_congruent_list_SCZ(:,i_SCZ)  , fixed_params_0);%, initial_params,0); 
        [param_fitted_incongruent,~,~]  = sigm_fit_weight_by_SD(relative_amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), SE_Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), fixed_params_0);%, initial_params,0);
                
        if abs(param_fitted_incongruent(3)>15) | param_fitted_incongruent(4)<0 | param_fitted_incongruent(4)>2.5
            [param_fitted_incongruent,~,~]  = sigm_fit(relative_amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ), fixed_params_0);%, initial_params,0);
        end
    end
    
    %%% Can generalize to MLE weighted by n_obs etc, but don't think that's needed.
    
    param_fit_sigm_list_SCZ(:,i_SCZ) = param_fitted;
    param_fit_congruent_sigm_list_SCZ(:,i_SCZ) = param_fitted_congruent;
    param_fit_incongruent_sigm_list_SCZ(:,i_SCZ) = param_fitted_incongruent;

    if flag_plot_SCZ==1
        %%% Plots
        figure;
        subplot(3,3,1)
        hold on;
        plot(net_evidence_smooth_list, sigmoid_fit_function_list(param_fitted));
        scatter(net_evidence_list, Pcorr_amb_list_SCZ(:,i_SCZ));
        ylabel('Proportion correct');
        xlabel('Net evidence');
        ylim([0.5,1]);
        title('All trials')

        subplot(3,3,4)
        hold on;
        plot(amb_smooth_list, sigmoid_fit_function_list(param_fitted));
        scatter(amb_list, Pcorr_amb_list_SCZ(:,i_SCZ));
        ylabel('Proportion correct');
        xlabel('Ambiguity');
        title('All trials')
        ylim([0.5,1]);

        subplot(3,3,7)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted));
        scatter(relative_amb_list, Pcorr_amb_list_SCZ(:,i_SCZ));
        ylabel('Proportion correct');
        xlabel('Relative ambiguity');
        title('All trials')
        ylim([0.5,1]);

        subplot(3,3,2)
        hold on;
        plot(net_evidence_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent));
        scatter(net_evidence_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ));
        ylabel('Proportion correct');
        xlabel('Net evidence');
        title('Incongruent trials')

        subplot(3,3,5)
        hold on;
        plot(amb_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent));
        scatter(amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ));
        ylabel('Proportion correct');
        xlabel('Ambiguity');
        title('Incongruent trials')

        subplot(3,3,8)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent));
        scatter(relative_amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ));
        ylabel('Proportion correct');
        xlabel('Relative ambiguity');
        title('Incongruent trials')

        subplot(3,3,3)
        hold on;
        plot(net_evidence_smooth_list, sigmoid_fit_function_list(param_fitted_congruent));
        scatter(net_evidence_list, Pcorr_amb_congruent_list_SCZ(:,i_SCZ));
        ylabel('Proportion correct');
        xlabel('Net evidence');
        title('Congruent trials')
        ylim([0.5,1]);

        subplot(3,3,6)
        hold on;
        plot(amb_smooth_list, sigmoid_fit_function_list(param_fitted_congruent));
        scatter(amb_list, Pcorr_amb_congruent_list_SCZ(:,i_SCZ));
        ylabel('Proportion correct');
        xlabel('Ambiguity');
        title('Congruent trials')
        ylim([0.5,1]);

        subplot(3,3,9)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted_congruent));
        scatter(relative_amb_list, Pcorr_amb_congruent_list_SCZ(:,i_SCZ));
        ylabel('Proportion correct');
        xlabel('Relative ambiguity');
        title('Congruent trials')
        ylim([0.5,1]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This code section visualizes the distribution of fitted parameters (ymax,
%  EC50, and slope) from the psychometric fits for both Healthy Controls
%  (HCS) and Schizophrenia patients (SCZ). It first selects subjects that
%  meet a criterion for correct performance on non-conflict trials,
%  then plots scatter or box plots comparing HCS vs. SCZ, as well as
%  subdividing Incongruent vs. Congruent conditions.
%
%  Key steps:
%  1. Subject selection based on performance criteria
%  2. Generation of scatter or boxplot diagrams for 'ymax', 'EC50', and 'slope'
%  3. Comparison between HCS and SCZ or between Congruent and Incongruent groups
%  4. Optional statistical tests (e.g., ranksum) to compare distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Distribution of fitted params

% Define a performance criterion for including subjects:
Pcorr_noconfl_criteria = 0.7; % temp for plotting performance overlayed .
Pcorr_noconfl_criteria_2 = 0.64;

% Identify which subjects meet the criterion for HCS and SCZ groups:
i_HCS_subs_to_use = find(Pcorr_amb_incongruent_list_HCS(1,:)>Pcorr_noconfl_criteria);
i_SCZ_subs_to_use = find(Pcorr_amb_incongruent_list_SCZ(1,:)>Pcorr_noconfl_criteria);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot distribution of fit parameters (ymax, EC50, slope)
%  for both HCS and SCZ, and for different trial types (All, Incongruent,
%  Congruent). The code uses arrays such as:
%    param_fit_sigm_list_HCS     : Fit parameters (all trials) for HCS
%    param_fit_incongruent_sigm_list_HCS : Fit parameters for incongruent, HCS
%    param_fit_congruent_sigm_list_HCS   : Fit parameters for congruent, HCS
%  Similarly for SCZ. The indexing (2,3,4) refers to param array elements:
%    2 = ymax, 3 = EC50, 4 = slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
%--------------------------------------------------------------------------
% Subplot (3,1,1) : Ymax for each group and condition
%--------------------------------------------------------------------------
subplot(3,1,1)
toplot{1} = param_fit_sigm_list_HCS(2,i_HCS_subs_to_use);
toplot{2} = param_fit_sigm_list_SCZ(2,i_SCZ_subs_to_use);
toplot{3} = param_fit_incongruent_sigm_list_HCS(2,i_HCS_subs_to_use);
toplot{4} = param_fit_incongruent_sigm_list_SCZ(2,i_SCZ_subs_to_use);
toplot{5} = param_fit_congruent_sigm_list_HCS(2,i_HCS_subs_to_use);
toplot{6} = param_fit_congruent_sigm_list_SCZ(2,i_SCZ_subs_to_use);
hold on 

% Scatter points for each grouping, at slightly different x-locations:
scatter(ones(length(toplot{1}),1), toplot{1})
scatter(4*ones(length(toplot{3}),1), toplot{3})
scatter(7*ones(length(toplot{5}),1), toplot{5})
scatter(2*ones(length(toplot{2}),1), toplot{2})
scatter(5*ones(length(toplot{4}),1), toplot{4})
scatter(8*ones(length(toplot{6}),1), toplot{6})
xticks([1,2,4,5,7,8])
xticklabels(["HCS", "SCZ", "HCS, incon", "SCZ, incon", "HCS, cong", "SCZ, cong"])
ylabel('ymax')
% ylim([0, 2])

%--------------------------------------------------------------------------
% Subplot (3,1,2) : EC50 for each group and condition
%--------------------------------------------------------------------------
subplot(3,1,2)
toplot{1} = param_fit_sigm_list_HCS(3,i_HCS_subs_to_use);
toplot{2} = param_fit_sigm_list_SCZ(3,i_SCZ_subs_to_use);
toplot{3} = param_fit_incongruent_sigm_list_HCS(3,i_HCS_subs_to_use);
toplot{4} = param_fit_incongruent_sigm_list_SCZ(3,i_SCZ_subs_to_use);
toplot{5} = param_fit_congruent_sigm_list_HCS(3,i_HCS_subs_to_use);
toplot{6} = param_fit_congruent_sigm_list_SCZ(3,i_SCZ_subs_to_use);
hold on 

scatter(ones(length(toplot{1}),1), toplot{1})
scatter(4*ones(length(toplot{3}),1), toplot{3})
scatter(7*ones(length(toplot{5}),1), toplot{5})
scatter(2*ones(length(toplot{2}),1), toplot{2})
scatter(5*ones(length(toplot{4}),1), toplot{4})
scatter(8*ones(length(toplot{6}),1), toplot{6})
xticks([1,2,4,5,7,8])
xticklabels(["HCS", "SCZ", "HCS, incon", "SCZ, incon", "HCS, cong", "SCZ, cong"])
ylabel('EC50')
% ylim([-1, 2])
ylim([-2, 6])

%--------------------------------------------------------------------------
% Subplot (3,1,3) : Slope for each group and condition
%--------------------------------------------------------------------------
subplot(3,1,3)
toplot{1} = param_fit_sigm_list_HCS(4,i_HCS_subs_to_use);
toplot{2} = param_fit_sigm_list_SCZ(4,i_SCZ_subs_to_use);
toplot{3} = param_fit_incongruent_sigm_list_HCS(4,i_HCS_subs_to_use);
toplot{4} = param_fit_incongruent_sigm_list_SCZ(4,i_SCZ_subs_to_use);
toplot{5} = param_fit_congruent_sigm_list_HCS(4,i_HCS_subs_to_use);
toplot{6} = param_fit_congruent_sigm_list_SCZ(4,i_SCZ_subs_to_use);
hold on 

scatter(ones(length(toplot{1}),1), toplot{1})
scatter(4*ones(length(toplot{3}),1), toplot{3})
scatter(7*ones(length(toplot{5}),1), toplot{5})
scatter(2*ones(length(toplot{2}),1), toplot{2})
scatter(5*ones(length(toplot{4}),1), toplot{4})
scatter(8*ones(length(toplot{6}),1), toplot{6})
xticks([1,2,4,5,7,8])
xticklabels(["HCS", "SCZ", "HCS, incon", "SCZ, incon", "HCS, cong", "SCZ, cong"])
ylabel('Slope')
% ylim([-7, 5])
ylim([0, 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next figure: Incongruent-trials-specific box/scatter plots comparing
% HCS vs. SCZ for each parameter. This figure uses 'boxplot2(...)' plus
% scatter overlays, focusing only on incongruent-trial fits. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
%--------------------------------------------------------------------------
% Subplot (3,2,1) : Boxplot for Ymax (incongruent) across HCS vs. SCZ
%--------------------------------------------------------------------------
subplot(3,2,1)
toplot{1} = param_fit_incongruent_sigm_list_HCS(2,i_HCS_subs_to_use);
toplot{2} = param_fit_incongruent_sigm_list_SCZ(2,i_SCZ_subs_to_use);
hold on 
boxplot2([1,2],toplot)%, 'BoxStyle','filled')
scatter(ones(length(toplot{1}),1), toplot{1})
scatter(2*ones(length(toplot{2}),1), toplot{2})
xticks([1,2])
xticklabels(["HCS", "SCZ"])
ylabel('ymax')
title('Incongruent trials')
% ylim([0, 2])
xlim([0.5, 2.5])

%--------------------------------------------------------------------------
% Subplot (3,2,3) : Boxplot for EC50 (incongruent) across HCS vs. SCZ
%--------------------------------------------------------------------------
subplot(3,2,3)
toplot{1} = param_fit_incongruent_sigm_list_HCS(3,i_HCS_subs_to_use);
toplot{2} = param_fit_incongruent_sigm_list_SCZ(3,i_SCZ_subs_to_use);
hold on 
boxplot2([1,2],toplot)%, 'BoxStyle','filled')
scatter(ones(length(toplot{1}),1), toplot{1})
scatter(2*ones(length(toplot{2}),1), toplot{2})
xticks([1,2])
xticklabels(["HCS", "SCZ"])
ylabel('EC50 [net evidence]')
title('Incongruent trials')
% ylim([-1, 2])
% ylim([-2, 6])
xlim([0.5, 2.5])

%--------------------------------------------------------------------------
% Subplot (3,2,5) : Boxplot for Slope (incongruent) across HCS vs. SCZ
%--------------------------------------------------------------------------
subplot(3,2,5)
toplot{1} = param_fit_incongruent_sigm_list_HCS(4,i_HCS_subs_to_use);
toplot{2} = param_fit_incongruent_sigm_list_SCZ(4,i_SCZ_subs_to_use);
hold on 
boxplot2([1,2],toplot)%, 'BoxStyle','filled')
scatter(ones(length(toplot{1}),1), toplot{1})
scatter(2*ones(length(toplot{2}),1), toplot{2})
xticks([1,2])
xticklabels(["HCS", "SCZ"])
ylabel('Slope [1/net evidence]')
title('Incongruent trials')
% ylim([-7, 5])
% ylim([0, 2])
xlim([0.5, 2.5])

%--------------------------------------------------------------------------
% Subplot (3,2,4) : Boxplot for EC50 in "relative ambiguity" terms
%--------------------------------------------------------------------------
subplot(3,2,4)
% We transform the net evidence param for each subject into an equivalent
% "relative ambiguity" measure. Then apply a rank-sum test:
toplot{1} = (n_non0_pulse - param_fit_incongruent_sigm_list_HCS(3,i_HCS_subs_to_use)) ./ (n_non0_pulse + param_fit_incongruent_sigm_list_HCS(3,i_HCS_subs_to_use));
toplot{2} = (n_non0_pulse - param_fit_incongruent_sigm_list_SCZ(3,i_SCZ_subs_to_use)) ./ (n_non0_pulse + param_fit_incongruent_sigm_list_SCZ(3,i_SCZ_subs_to_use));
% Clean out invalid values (<0):
EC50_amb_HCS = (n_non0_pulse - param_fit_incongruent_sigm_list_HCS(3,i_HCS_subs_to_use)) ./ (n_non0_pulse + param_fit_incongruent_sigm_list_HCS(3,i_HCS_subs_to_use));
EC50_amb_SCZ = (n_non0_pulse - param_fit_incongruent_sigm_list_SCZ(3,i_SCZ_subs_to_use)) ./ (n_non0_pulse + param_fit_incongruent_sigm_list_SCZ(3,i_SCZ_subs_to_use));

EC50_amb_HCS(EC50_amb_HCS<0) = NaN;
toplot{1}(toplot{1}<0) = NaN;
EC50_amb_SCZ(EC50_amb_SCZ<0) = NaN;
toplot{2}(toplot{2}<0) = NaN;

% Conduct a rank-sum test comparing the two groups:
p_ranksum_EC50_HCS_SCZ = ranksum(EC50_amb_HCS, EC50_amb_SCZ);

hold on 
boxplot2([1,2],toplot)%, 'BoxStyle','filled')
scatter(ones(length(toplot{1}),1)+0.1*randn(length(toplot{1}),1), toplot{1})
scatter(2*ones(length(toplot{2}),1)+0.1*randn(length(toplot{2}),1), toplot{2})
xticks([1,2])
xticklabels(["HCS", "SCZ"])
ylabel('EC50 [relative amb]')
title('Incongruent trials')
% ylim([-1, 2])
% ylim([-2, 6])
xlim([0.5, 2.5])


%--------------------------------------------------------------------------
% Subplot (3,2,6) : Another measure of slope in relative-amb terms
%--------------------------------------------------------------------------
subplot(3,2,6)
%%% N.B.: need to verify. Compare value to when I fit over rel amb (but with modified sigmoid (1/(1+10^(1/x))) etc).
toplot{1} = 1./((n_non0_pulse - 1./param_fit_incongruent_sigm_list_HCS(4,i_HCS_subs_to_use)) ./ (n_non0_pulse + 1./param_fit_incongruent_sigm_list_HCS(4,i_HCS_subs_to_use)));
toplot{2} = 1./((n_non0_pulse - 1./param_fit_incongruent_sigm_list_SCZ(4,i_SCZ_subs_to_use)) ./ (n_non0_pulse + 1./param_fit_incongruent_sigm_list_SCZ(4,i_SCZ_subs_to_use)));
hold on 
boxplot2([1,2],toplot)%, 'BoxStyle','filled')
scatter(ones(length(toplot{1}),1)+0.1*randn(length(toplot{1}),1), toplot{1})
scatter(2*ones(length(toplot{2}),1)+0.1*randn(length(toplot{2}),1), toplot{2})
xticks([1,2])
xticklabels(["HCS", "SCZ"])
ylabel('Slope [1/relative amb]')
title('Incongruent trials')
% ylim([-7, 5])
% ylim([0, 2])
xlim([0.5, 2.5])
p_ranksum_slope_HCS_SCZ = ranksum(toplot{1}, toplot{2});

% Save figure with the specified renderer to PNG/PDF/FIG:
set(gcf,'Renderer','painters')
saveas(gcf, strcat('fit_param_box_scatter_incongruent.png'));
saveas(gcf, strcat('fit_param_box_scatter_incongruent.pdf'));
saveas(gcf, strcat('fit_param_box_scatter_incongruent.fig'));

%% Plot Performance.

if flag_fitting_over_rel_amb==0
    [param_fitted_incongruent_avg_HCS,stats_incongruent_avg_HCS,~]  = sigm_fit(repmat(net_evidence_list',length(i_HCS_subs_to_use),1), reshape(Pcorr_amb_incongruent_list_HCS(:,i_HCS_subs_to_use),[],1), fixed_params_0);%, initial_params,0);     %%%%% Fit over subjects' performance.
    %%%%% could also fit over each trial, but would need loop above to output response_list... not worth.
elseif flag_fitting_over_rel_amb==1
    %%% Fit over relative ambiguity instead
    [param_fitted_incongruent_avg_HCS,stats_incongruent_avg_HCS,~]  = sigm_fit(repmat(relative_amb_list',length(i_HCS_subs_to_use),1), reshape(Pcorr_amb_incongruent_list_HCS(:,i_HCS_subs_to_use),[],1), fixed_params_0);%, initial_params,0);     %%%%% Fit over subjects' performance.
    %%%%% could also fit over each trial, but would need loop above to output response_list... not worth.
end
% response_list_incongruent

% --- std & CI calculation for HCS incongruent performance ---
% We compute the standard deviation across subjects for each row 
% (i.e., each ambiguity level), then get a T-based confidence interval 
% before ultimately converting it to SEM.
std_Pcorr_amb_incongruent_HCS = std(Pcorr_amb_incongruent_list_HCS(:,i_HCS_subs_to_use),[],2);
CI_Pcorr_amb_incongruent_HCS = tinv([0.975],length(i_HCS_subs_to_use)-1)'.*(std_Pcorr_amb_incongruent_HCS)/length(i_HCS_subs_to_use)^0.5;
% Convert that CI to SEM:
CI_Pcorr_amb_incongruent_HCS = (std_Pcorr_amb_incongruent_HCS)/length(i_HCS_subs_to_use)^0.5;

figure;

% -------------------
% 1) HCS group, incongruent trials
%    Subplot (2,2,1) : Fitted curve + raw means and SEM
subplot(2,2,1)
hold on;
plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent_avg_HCS));
shadedErrorBar(relative_amb_list, mean(Pcorr_amb_incongruent_list_HCS(:,i_HCS_subs_to_use),2) , CI_Pcorr_amb_incongruent_HCS, 'lineProps',{'linewidth',1.5})
ylim([0.5, 1])
ylabel('Proportion correct');
xlabel('Relative ambiguity');
title('HCS, Incongruent trials')

% -------------------
% 2) HCS group, incongruent trials
%    Subplot (2,2,2): Display curve with confidence intervals from stats_incongruent_avg_HCS
subplot(2,2,2)
hold on;
shadedErrorBar(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent_avg_HCS), (stats_incongruent_avg_HCS.ypred - stats_incongruent_avg_HCS.ypredlowerCI)/2, 'lineProps',{'linewidth',1.5,'color','k'})
% Also scatter the raw means for reference
scatter(relative_amb_list, mean(Pcorr_amb_incongruent_list_HCS(:,i_HCS_subs_to_use),2),'k');
ylim([0.5, 1])
ylabel('Proportion correct');
xlabel('Relative ambiguity');
title('HCS, Incongruent trials')

% --- Now do the same for SCZ group, aggregated fits ---
if flag_fitting_over_rel_amb==0
    % Fitting SCZ data using net evidence:
    [param_fitted_incongruent_avg_SCZ,stats_incongruent_avg_SCZ,~]  = sigm_fit(repmat(net_evidence_list',length(i_SCZ_subs_to_use),1), reshape(Pcorr_amb_incongruent_list_SCZ(:,i_SCZ_subs_to_use),[],1), fixed_params_0);%, initial_params,0);     %%%%% Fit over subjects' performance.
    %%%%% could also fit over each trial, but would need loop above to output response_list... not worth.
elseif flag_fitting_over_rel_amb==1
    % Fitting SCZ data using relative ambiguity:
    [param_fitted_incongruent_avg_SCZ,stats_incongruent_avg_SCZ,~]  = sigm_fit(repmat(relative_amb_list',length(i_SCZ_subs_to_use),1), reshape(Pcorr_amb_incongruent_list_SCZ(:,i_SCZ_subs_to_use),[],1), fixed_params_0);%, initial_params,0);     %%%%% Fit over subjects' performance.
    %%%%% could also fit over each trial, but would need loop above to output response_list... not worth.
end

% --- std & CI for SCZ incongruent performance ---
std_Pcorr_amb_incongruent_SCZ = std(Pcorr_amb_incongruent_list_SCZ(:,i_SCZ_subs_to_use),[],2);
CI_Pcorr_amb_incongruent_SCZ = tinv([0.975],length(i_SCZ_subs_to_use)-1)'.*(std_Pcorr_amb_incongruent_SCZ)/length(i_SCZ_subs_to_use)^0.5;       
% Convert that CI to SEM:
CI_Pcorr_amb_incongruent_SCZ = (std_Pcorr_amb_incongruent_SCZ)/length(i_SCZ_subs_to_use)^0.5;       

% -------------------
% 3) SCZ group, incongruent trials
%    Subplot (2,2,3): SCZ curve + raw means + SEM
subplot(2,2,3)
hold on;
plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent_avg_SCZ));
shadedErrorBar(relative_amb_list, mean(Pcorr_amb_incongruent_list_SCZ(:,i_SCZ_subs_to_use),2) , CI_Pcorr_amb_incongruent_SCZ, 'lineProps',{'linewidth',1.5})
ylim([0.5, 1])
ylabel('Proportion correct');
xlabel('Relative ambiguity');
title('SCZ, Incongruent trials')

% -------------------
% 4) SCZ group, incongruent trials
%    Subplot (2,2,4): SCZ aggregated fit with confidence intervals
subplot(2,2,4)
hold on;
shadedErrorBar(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent_avg_SCZ), (stats_incongruent_avg_SCZ.ypred - stats_incongruent_avg_SCZ.ypredlowerCI)/2, 'lineProps',{'linewidth',1.5,'color',[227,26,28]/255})

scatter(relative_amb_list, mean(Pcorr_amb_incongruent_list_SCZ(:,i_SCZ_subs_to_use),2));

%%% Overlap SCZ plot with HCS fit for direct comparison:
shadedErrorBar(relative_amb_smooth_list, sigmoid_fit_function_list(param_fitted_incongruent_avg_HCS), (stats_incongruent_avg_HCS.ypred - stats_incongruent_avg_HCS.ypredlowerCI)/2, 'lineProps',{'linewidth',1.5,'color','k'})
%%%%% temp
scatter(relative_amb_list, mean(Pcorr_amb_incongruent_list_HCS(:,i_HCS_subs_to_use),2),'k');

ylim([0.5, 1])
ylabel('Proportion correct');
xlabel('Relative ambiguity');
title('SCZ, Incongruent trials')

% Save figure in multiple formats
set(gcf,'Renderer','painters')
saveas(gcf, strcat('Pcorr_over_subjects.png'));
saveas(gcf, strcat('Pcorr_over_subjects.pdf'));
saveas(gcf, strcat('Pcorr_over_subjects.fig'));


%% Stat test
% Here we compare HCS vs. SCZ performance at each ambiguity level:
%   1) ranksum test: checks whether the medians differ between HCS vs. SCZ 
%      (non-parametric).
%   2) chi2test: compares success/failure frequencies across groups 
%      at that ambiguity level.
p_ranksum_HCS_SCZ = zeros(4,1);
p_chi2_HCS_SCZ = zeros(4,1);

% n_trials_incongruent_each / n_trials_congruent_each: 
% approximate number of trials used to form success/failure
n_trials_incongruent_each = 80;
n_trials_congruent_each = 40;

for i_ne = 1:4
    % ranksum compares distributions of Pcorr_amb_incongruent_list_HCS vs. SCZ
    p_ranksum_HCS_SCZ(i_ne) = ranksum(Pcorr_amb_incongruent_list_HCS(i_ne,i_HCS_subs_to_use), Pcorr_amb_incongruent_list_SCZ(i_ne,i_SCZ_subs_to_use));       
    % chi2test tests a 2×2 contingency table for success/fail in HCS vs. SCZ
    % success/fail counts are approximated by Pcorr * #trials vs. (1 - Pcorr)*#trials
    [p_chi2_HCS_SCZ(i_ne),~] = chi2test([sum(n_trials_incongruent_each * Pcorr_amb_incongruent_list_HCS(i_ne,i_HCS_subs_to_use)), sum(n_trials_incongruent_each * (1-Pcorr_amb_incongruent_list_HCS(i_ne,i_HCS_subs_to_use)))  ;  sum(n_trials_incongruent_each * Pcorr_amb_incongruent_list_SCZ(i_ne,i_SCZ_subs_to_use)), sum(n_trials_incongruent_each * (1-Pcorr_amb_incongruent_list_SCZ(i_ne,i_SCZ_subs_to_use)))]);

end

%% Plot HCS and SCZ performance plots, overlayed.
% -------------------------------------------------------------------------
% In this section, we create two figures (101 for HCS, 102 for SCZ) showing 
% individual subject curves overlaid in subplots. Each subject’s data is 
% plotted as a line + points for raw performance. We also distinguish 
% "type 1" vs. "type 2" subjects based on how much their performance drops 
% from no conflict to moderate conflict trials.
% -------------------------------------------------------------------------
type_1_2_thres_temp = 0.1;
ylim_temp = [0.5,1];

% Check if we want to generate overlay plots:
if flag_plot_overlayed==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure 101: Overlays for HCS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(101);
    for ind_HCS = 1:length(i_HCS_subs_to_use) % if use neurons that fit criteria (generally, performance at lowest conflict).
        % i_HCS_subs_to_use: indices of HCS subjects who met a performance 
        % criterion on low-conflict trials.
        i_HCS = i_HCS_subs_to_use(ind_HCS);


        %------------------------------------------------------------------
        % Subplot (2,3,1) - raw incongruent curve for HCS, 
        % plotting proportion correct vs. relative ambiguity.
        %------------------------------------------------------------------
        subplot(2,3,1)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_HCS(:,i_HCS)), 'Color',[1,1,1]*0.5*(1+i_HCS/length(data_list_HCS)) );
        scatter(relative_amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS));
        ylim(ylim_temp)
        ylabel('Proportion correct');
        xlabel('Relative ambiguity');
        title('Incongruent trials')
            
        %------------------------------------------------------------------
        % Subplot (2,3,4) - same as above, but normalized to that 
        % subject's maximum performance (i.e., dividing by the 
        % subject’s performance at the easiest condition).
        % This highlights relative performance changes across ambiguity.
        %------------------------------------------------------------------
        subplot(2,3,4)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_HCS(:,i_HCS))/Pcorr_amb_incongruent_list_HCS(1,i_HCS), 'Color',[1,1,1]*0.5*(1+i_HCS/length(data_list_HCS)) );
        scatter(relative_amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS)/Pcorr_amb_incongruent_list_HCS(1,i_HCS));
        ylabel('Proportion correct, relative to asymptote');
        xlabel('Relative ambiguity');
        title('Incongruent trials')            

        %------------------------------------------------------------------
        % Distinguish "type 1" vs. "type 2" subjects based on whether the 
        % drop from no conflict (lowest conflict) to moderate conflict 
        % is less than a threshold (type_1_2_thres_temp).
        % Place them in separate subplots for further visualization.
        %------------------------------------------------------------------
        if nanmax(Pcorr_amb_incongruent_list_HCS(1:2,i_HCS))-Pcorr_amb_incongruent_list_HCS(3,i_HCS)<type_1_2_thres_temp
            % Subplot (2,3,2) & (2,3,5) for "type 1" subjects
            subplot(2,3,2)
            hold on;
            plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_HCS(:,i_HCS)), 'Color',[1,1,1]*0.5*(1+i_HCS/length(data_list_HCS)) );
            scatter(relative_amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS));
            ylim(ylim_temp)
            ylabel('Proportion correct');
            xlabel('Relative ambiguity');
            title('Incongruent trials, type 1 subjects')            
            
            subplot(2,3,5)
            hold on;
            plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_HCS(:,i_HCS))/Pcorr_amb_incongruent_list_HCS(1,i_HCS), 'Color',[1,1,1]*0.5*(1+i_HCS/length(data_list_HCS)) );
            scatter(relative_amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS)/Pcorr_amb_incongruent_list_HCS(1,i_HCS));
            ylabel('Proportion correct, relative to asymptote');
            xlabel('Relative ambiguity');
            title('Incongruent trials, type 1 subjects')            
        else
            % Subplot (2,3,3) & (2,3,6) for "type 2" subjects
            subplot(2,3,3)
            hold on;
            plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_HCS(:,i_HCS)), 'Color',[1,1,1]*0.5*(1+i_HCS/length(data_list_HCS)) );
            scatter(relative_amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS));
            ylim(ylim_temp)
            ylabel('Proportion correct');
            xlabel('Relative ambiguity');
            title('Incongruent trials, type 2 subjects')        
            
            subplot(2,3,6)
            hold on;
            plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_HCS(:,i_HCS))/Pcorr_amb_incongruent_list_HCS(1,i_HCS), 'Color',[1,1,1]*0.5*(1+i_HCS/length(data_list_HCS)) );
            scatter(relative_amb_list, Pcorr_amb_incongruent_list_HCS(:,i_HCS)/Pcorr_amb_incongruent_list_HCS(1,i_HCS));
            ylabel('Proportion correct, relative to asymptote');
            xlabel('Relative ambiguity');
            title('Incongruent trials, type 2 subjects')                        
        end
        
    end
    sgtitle('HCS') % Overall title for the figure
    set(gcf,'Renderer','painters')
    set(gcf,'Position',get(0,'Screensize')); % Maximize figure window
    
    % Save the HCS overlay figure in multiple formats
    saveas(gcf, strcat('Pcorr_subjects_overlayed_HCS.png'));
    saveas(gcf, strcat('Pcorr_subjects_overlayed_HCS.pdf'));
    saveas(gcf, strcat('Pcorr_subjects_overlayed_HCS.fig'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure 102: Overlays for SCZ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop over SCZ subjects who passed a performance criterion on 
    % low-conflict trials.
    figure(102);
    for ind_SCZ = 1:length(i_SCZ_subs_to_use) % if use neurons that fit criteria (generally, performance at lowest conflict).
        i_SCZ = i_SCZ_subs_to_use(ind_SCZ);
            
        %----------------------------------------------------------------------
        % Subplot (2,3,1): Plot raw incongruent performance for each SCZ subject
        % vs. relative ambiguity. We also show the subject-specific fitted curve.
        %----------------------------------------------------------------------
        subplot(2,3,1)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_SCZ(:,i_SCZ)), 'Color',[1,1,1]*0.5*(1+i_SCZ/length(data_list_SCZ)) );
        % Scatter raw data points:
        scatter(relative_amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ));
        ylim(ylim_temp)
        ylabel('Proportion correct');
        xlabel('Relative ambiguity');
        title('Incongruent trials')
        
        %----------------------------------------------------------------------
        % Subplot (2,3,4): Normalized performance for SCZ 
        % (dividing by the subject’s no-conflict or easiest condition).
        %----------------------------------------------------------------------        
        subplot(2,3,4)
        hold on;
        plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_SCZ(:,i_SCZ))/Pcorr_amb_incongruent_list_SCZ(1,i_SCZ), 'Color',[1,1,1]*0.5*(1+i_SCZ/length(data_list_SCZ)) );
        % Scatter raw points, also normalized by easiest condition:
        scatter(relative_amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ)/Pcorr_amb_incongruent_list_SCZ(1,i_SCZ));
        ylabel('Proportion correct, relative to asymptote');
        xlabel('Relative ambiguity');
        title('Incongruent trials')            

        %------------------------------------------------------------------
        % Classify SCZ subject as "type 1" or "type 2", same logic as HCS
        %------------------------------------------------------------------
        if nanmax(Pcorr_amb_incongruent_list_SCZ(1:2,i_SCZ))-Pcorr_amb_incongruent_list_SCZ(3,i_SCZ)<type_1_2_thres_temp
            % Subplot (2,3,2) & (2,3,5) for "type 1" SCZ
            subplot(2,3,2)
            hold on;
            plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_SCZ(:,i_SCZ)), 'Color',[1,1,1]*0.5*(1+i_SCZ/length(data_list_SCZ)) );
            scatter(relative_amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ));
            ylim(ylim_temp)
            ylabel('Proportion correct');
            xlabel('Relative ambiguity');
            title('Incongruent trials, type 1 subjects')            
            
            subplot(2,3,5)
            hold on;
            plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_SCZ(:,i_SCZ))/Pcorr_amb_incongruent_list_SCZ(1,i_SCZ), 'Color',[1,1,1]*0.5*(1+i_SCZ/length(data_list_SCZ)) );
            scatter(relative_amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ)/Pcorr_amb_incongruent_list_SCZ(1,i_SCZ));
            ylabel('Proportion correct, relative to asymptote');
            xlabel('Relative ambiguity');
            title('Incongruent trials, type 1 subjects')            

        else
            % Subplot (2,3,3) & (2,3,6) for "type 2" SCZ
            subplot(2,3,3)
            hold on;
            plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_SCZ(:,i_SCZ)), 'Color',[1,1,1]*0.5*(1+i_SCZ/length(data_list_SCZ)) );
            scatter(relative_amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ));
            ylim(ylim_temp)
            ylabel('Proportion correct');
            xlabel('Relative ambiguity');
            title('Incongruent trials, type 2 subjects')        
            
            subplot(2,3,6)
            hold on;
            plot(relative_amb_smooth_list, sigmoid_fit_function_list(param_fit_incongruent_sigm_list_SCZ(:,i_SCZ))/Pcorr_amb_incongruent_list_SCZ(1,i_SCZ), 'Color',[1,1,1]*0.5*(1+i_SCZ/length(data_list_SCZ)) );
            scatter(relative_amb_list, Pcorr_amb_incongruent_list_SCZ(:,i_SCZ)/Pcorr_amb_incongruent_list_SCZ(1,i_SCZ));
            ylabel('Proportion correct, relative to asymptote');
            xlabel('Relative ambiguity');
            title('Incongruent trials, type 2 subjects')                        
        end
        
    end
    % Annotate the entire figure for SCZ
    sgtitle('SCZ')
    set(gcf,'Renderer','painters')
    set(gcf,'Position',get(0,'Screensize'));
    
    % Save figure for SCZ overlays  
    saveas(gcf, strcat('Pcorr_subjects_overlayed_SCZ.png'));
    saveas(gcf, strcat('Pcorr_subjects_overlayed_SCZ.pdf'));
    saveas(gcf, strcat('Pcorr_subjects_overlayed_SCZ.fig'));
end

%% Example for model goodness-of-fit. Log-likelihood etc
%%% Several assumptions here... I'm only considering subjects used in the code above (i_HCS_subs_to_use, i_SCZ_subs_to_use).
% Also assuming for simplicity that n_total_trials are the same across all conditions of net evidence levels. If not. the generalization is LL = sum_i=1,3,5,9 n_trials_ne=i * sum(p_data_ne=i * log(p_model_ne=i) + (1-p_data_ne=i) * log(1-p_model_ne=i)) and repeat over subjects.
%%% Finally, note that the values are meaningless unless compared to that from other models.

LL_HCS_sigm_fit = 0;
LL_SCZ_sigm_fit = 0;

% Accumulate log-likelihood for HCS, using their param_fitted_incongruent_avg_HCS
for i_temp = 1:length(i_HCS_subs_to_use)
    LL_HCS_sigm_fit = LL_HCS_sigm_fit + sum(Pcorr_amb_incongruent_list_HCS(:,i_HCS_subs_to_use(i_temp)) .* log(sigmoid_fit_function_list_2(param_fitted_incongruent_avg_HCS))') + sum((1-Pcorr_amb_incongruent_list_HCS(:,i_HCS_subs_to_use(i_temp))) .* log(1-sigmoid_fit_function_list_2(param_fitted_incongruent_avg_HCS))');
end

% AIC/BIC: we have 3 free parameters in the sigmoidal model 
% (since y-min was fixed at 0.5, the rest are free).
[aic_HCS_sigm_fit, bic_HCS_sigm_fit] = aicbic(LL_HCS_sigm_fit, 3, n_trials_task*length(i_HCS_subs_to_use)); % 2nd input argument = number of 3 params. e.g. = 3 if using sigm_fit with fixed params are [0.5, NaN, NaN, NaN].


% Accumulate log-likelihood for SCZ, using param_fitted_incongruent_avg_SCZ
for i_temp = 1:length(i_SCZ_subs_to_use)
    LL_SCZ_sigm_fit = LL_SCZ_sigm_fit + sum(Pcorr_amb_incongruent_list_SCZ(:,i_SCZ_subs_to_use(i_temp)) .* log(sigmoid_fit_function_list_2(param_fitted_incongruent_avg_SCZ))') + sum((1-Pcorr_amb_incongruent_list_SCZ(:,i_SCZ_subs_to_use(i_temp))) .* log(1-sigmoid_fit_function_list_2(param_fitted_incongruent_avg_SCZ))');
end
[aic_SCZ_sigm_fit, bic_SCZ_sigm_fit] = aicbic(LL_SCZ_sigm_fit, 3, n_trials_task*length(i_SCZ_subs_to_use)); % 2nd input argument = number of 3 params. e.g. = 3 if using sigm_fit with fixed params are [0.5, NaN, NaN, NaN].

%% Proportion of trials excluded per subject.
% This section plots how many trials were excluded based on each subject’s
% performance bouts and conflict conditions, highlighting differences 
% between HCS and SCZ groups.

figure;
subplot(1,1,1)
toplot{1} = P_trials_reject_list_HCS(i_HCS_subs_to_use);
toplot{2} = P_trials_reject_list_SCZ(i_SCZ_subs_to_use);
hold on 
boxplot2([1,2],toplot)%, 'BoxStyle','filled')
scatter(ones(length(toplot{1}),1), toplot{1})
scatter(2*ones(length(toplot{2}),1), toplot{2})
xticks([1,2])
xticklabels(["HCS", "SCZ"])
ylabel('Fraction trials excluded')
title('Incongruent trials')
% ylim([0, 2])
xlim([0.5, 2.5])

% Below, we collect easy/hard trial rejections for potential downstream use
toplot{3} = P_easy_trials_reject_list_HCS(i_HCS_subs_to_use);
toplot{4} = P_easy_trials_reject_list_SCZ(i_SCZ_subs_to_use);
toplot{5} = P_hard_trials_reject_list_HCS(i_HCS_subs_to_use);
toplot{6} = P_hard_trials_reject_list_SCZ(i_SCZ_subs_to_use);

% Save figure to various formats (PNG, PDF, and FIG)
set(gcf,'Renderer','painters')
saveas(gcf, strcat('P_trials_excluded_incongrueent.png'));
saveas(gcf, strcat('P_trials_excluded_incongrueent.pdf'));
saveas(gcf, strcat('P_trials_excluded_incongrueent.fig'));
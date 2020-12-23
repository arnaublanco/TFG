function [input_file_info] = getFileInfo(subject, Hem, CondClass, POIfile_ind)

% INPUT:
    % subject: Subject number
    % Hem: 'LH' or 'RH'
    % CondClass: Conditions to classify
    % POIfile_ind: POI file that will be used (1, 2 or 3)
    
% OUTPUT:
    % input_file_info: MATLAB object containing the data info.

dir_name = 'C:/Users/Arnau/Desktop/TFG/output/'; % Path to the files
dir_rois = 'C:/Users/Arnau/Desktop/TFG/ROIs/'; % Path to the ROIs/POIs
    
%%% Important parameters
nVols = 218;
nPreds = 18+1;  %% number of stimulation blocks plus 1 for baseline 
%CondClass=1:3; %% the conditions to classify, their codes in txt file
nTrials = 18;   %number of stimulation blocks (without baseline) per run
nPerRun = 6;    %% number of blocks per condition per run
nRuns = 4;

% !! The first subject underwent a short version of the experiment 

func_name = cell(1,nRuns);
dm_block = cell(1,nRuns);
dm_cond = cell(1,nRuns);
cond_locs_name = cell(1,4);

% !! Hem not used yet - i'm not separating by left and right hemisphere

for run = 1:nRuns
    
    if subject == 1
        w = 'short';
    else
        w = 'long';
    end
    
    s = "sub-0" + string(subject);
    file = s + "_ses-mri_task-AVScenes" + w +"_run-" + string(run) + "_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz";
    
    func_name{run} = dir_name + s + "/" + file; % Functional image file
    
    % File that maps single trials or blocks onto experimental conditions.
    cond_locs_name{run} = dir_name + s + "/" + s + "_ses-mri_task-AVScenes" + w + "_run-0" + run + "_events.tsv"; % EL_run1_AVScenesBlind_trialseq.txt
    
    % Directories for the design matrices (single block-wise and
    % condition-wise)
    dir_dm_b = dir_name + s + "/" + s + "_ses-mri_task-AVScenes" + w + "_run-0" + run + '_block.txt';
    dir_dm_c = dir_name + s + "/" + s + "_ses-mri_task-AVScenes" + w + "_run-0" + run + '_cond.txt';
    
    % If the design matrix files do not exist, create them in createDesignMatrix.m
    if ~exist(dir_dm_b,'file') || ~exist(dir_dm_c,'file')
        createDesignMatrix(subject,nPreds,cond_locs_name{run},dir_dm_b,dir_dm_c);
    end
    
    dm_block{run} = dir_dm_b;       % Design matrix (single trial or single block wise)
    dm_cond{run} = dir_dm_c;        % Design matrix (condition wise)
    
    % POI: Patch-of-Interest
    if POIfile_ind == 1
        poi_name = dir_rois + 'Auditory.nii'; % Auditory.poi
    elseif POIfile_ind == 2
        poi_name = dir_rois + 'Motor.nii'; % Motor.poi
    else
        poi_name = dir_rois + 'Visual.nii'; % Visual.poi
    end
    
end

%%% MATLAB object definition
input_file_info.dir_name = dir_name;
input_file_info.poi_name = poi_name;
input_file_info.func_name = func_name;
input_file_info.dm_name = dm_block;    % Block or trial wise
input_file_info.dm2 = dm_cond;         % Condition wise
input_file_info.CondClass = CondClass;
input_file_info.subject = subject;
input_file_info.cond_locs = cond_locs_name;
%input_file_info.Hem = Hem;

pars(1) = nVols;
pars(2) = nPreds;  %% remember to add 1 for the constant column
pars(3) = nTrials;
pars(4) = nPerRun;

input_file_info.pars = pars;
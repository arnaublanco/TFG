% Function that returns returns data information of the asked file.
%  INPUT:
%   · subject: Subject number
%   · Hem: 'LH' or 'RH'
%   · CondClass: Conditions to classify
%   · POIfile_ind: Index in POI file (1 -> Auditory, 2 -> Motor or 3 -> EVC)
%  OUTPUT:
%   · input_file_info: MATLAB object containing the data info.

function [input_file_info] = getFileInfo(subject, Hem, CondClass, POIfile_ind)

dir_name = '/Users/blancoarnau/Documents/ds002715/derivatives/fmriprep/'; % Path to the files
dir_rois = '/Users/blancoarnau/Documents/ds002715/derivatives/freesurfer/'; % Path to the ROIs/POIs
dir_events = '/Users/blancoarnau/Documents/ds002715/'; % Path to events.tsv

% Less volumes are selected because there's a glitch at the beginning of the signal.
if subject == 1
    nVols = 113; % The first subject underwent a shorter version of the experiment 
else
    nVols = 218;
end

nPreds = 18+1;  % Number of stimulation blocks (plus 1 for baseline) 
%CondClass=1:3; % Conditions to classify, their codes in txt file
nTrials = 18;   % Number of stimulation blocks (without baseline) per run
nPerRun = 6;    % Number of blocks per condition per run
nRuns = 4;      % Runs per subject

func_name = cell(1,nRuns);
dm_block = cell(1,nRuns);
dm_cond = cell(1,nRuns);
cond_locs_name = cell(1,4);

if strcmp(Hem,'LH')
    h = 'L';
else
    h = 'R';
end

if subject == 1
    w = 'short';
else
    w = 'long';
end

s = 'sub-0' + string(subject);

anat_name = dir_name + s + '/ses-mri/anat/' + s + '_ses-mri_run-1_hemi-' + h + '_inflated.surf.gii';

for run = 1:nRuns
    
    file = s + '_ses-mri_task-AVScenes' + w + '_run-' + string(run) + '_space-fsnative_hemi-' + h + '_bold.func.gii';
    
    func_name{run} = dir_name + s + '/ses-mri/func/' + file; % Functional image file
    
    % File that maps single trials or blocks onto experimental conditions.
    cond_locs_name{run} = dir_events + s + '/ses-mri/func/' + s + '_ses-mri_task-AVScenes' + w + '_run-0' + run + '_events.tsv'; % EL_run1_AVScenesBlind_trialseq.txt
    
    % Directories for the design matrices (single block-wise and condition-wise)
    dir_dm_b = dir_name + s + '/ses-mri/' + s + '_ses-mri_task-AVScenes' + w + '_run-0' + run + '_block.txt';
    dir_dm_c = dir_name + s + '/ses-mri/' + s + '_ses-mri_task-AVScenes' + w + '_run-0' + run + '_cond.txt';
    
    % If the design matrix files do not exist, create them in createDesignMatrix.m
    if ~exist(dir_dm_b,'file') || ~exist(dir_dm_c,'file')
        createDesignMatrix(subject,nPreds,cond_locs_name{run},dir_dm_b,dir_dm_c);
    end
    
    dm_block{run} = dir_dm_b;       % Design matrix (single trial or single block wise)
    dm_cond{run} = dir_dm_c;        % Design matrix (condition wise)
    
    % POI: Patch-of-Interest (or ROI: Region-of-Interest)
    if POIfile_ind == 1 % (to change laterrr!!!!!)
        poi_name = dir_rois + s + '/label/' + lower(Hem) + '.BA1_exvivo.label'; % Auditory (I'M GONNA USE BA1 BUT I HAVE TO CHANGE THIS LATER!!)
    elseif POIfile_ind == 2
        poi_name = dir_rois + s + '/label/' + lower(Hem) + '.BA4a_exvivo.label'; % Motor (primary motor area)
    else
        poi_name = [{dir_rois + s + '/label/' + lower(Hem) + '.V1_exvivo.label'},{dir_rois + s + '/label/' + lower(Hem) + '.V2_exvivo.label'}]; % Visual (V1 and V2)
    end
    
end

% MATLAB object definition
input_file_info.dir_name = dir_name + s + '/';
input_file_info.poi_name = poi_name;
input_file_info.func_name = func_name;
input_file_info.dm_name = dm_block;    % Block-wise
input_file_info.dm2 = dm_cond;         % Condition-wise
input_file_info.CondClass = CondClass;
input_file_info.subject = subject;
input_file_info.cond_locs = cond_locs_name;
input_file_info.Hem = Hem;
input_file_info.anat_name = anat_name;

pars(1) = nVols;
pars(2) = nPreds;  % Remember to add 1 for the constant column
pars(3) = nTrials;
pars(4) = nPerRun;

input_file_info.pars = pars;
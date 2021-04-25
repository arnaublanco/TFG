function [input_file_info] = getFileInfo(subject, Hem, CondClass, POIfile_ind)

% INPUT:
    % subject: Subject number
    % Hem: 'LH' or 'RH'
    % CondClass: Conditions to classify
    % POIfile_ind: Index in POI file (1 -> Auditory, 2 -> Motor or 3 -> EVC)
    
% OUTPUT:
    % input_file_info: MATLAB object containing the data info.

dir_name = '/Users/blancoarnau/Google Drive/TFG Arnau Blanco/data/output/'; % Path to the files
dir_rois = 'Extra/'; % Path to the ROIs
    
% Fewer volumes are selected because there is a glitch at the beginning of the signal.
if subject == 1
   nVols = 113;
else
    nVols = 218;
end
nPreds = 18+1;  % Number of stimulation blocks plus 1 for baseline 
nTrials = 18;   % Number of stimulation blocks (without baseline) per run
nPerRun = 6;    % Number of blocks per condition per run
nRuns = 4;

func_name = cell(1,nRuns);
dm_block = cell(1,nRuns);
dm_cond = cell(1,nRuns);
cond_locs_name = cell(1,4);

if subject == 1
    w = 'short';
else
    w = 'long';
end

for run = 1:nRuns
    
    % Functional image file
    func_name{run} = [dir_name,'sub-0',int2str(subject),'/sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-',int2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'];
    
    % File that maps single trials or blocks onto experimental conditions.
    cond_locs_name{run} = [dir_name,'sub-0',int2str(subject),'/sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-0',int2str(run),'_events.tsv'];
    
    % Directories for the design matrices (single block-wise and condition-wise)
    dir_dm_b = [dir_name,'sub-0',int2str(subject),'/sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-0',int2str(run),'_block.xls'];
    dir_dm_c = [dir_name,'sub-0',int2str(subject),'/sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-0',int2str(run),'_cond.xls'];
    
    % If the design matrix files do not exist, create them in createDesignMatrix.m
    if ~exist(dir_dm_b,'file') || ~exist(dir_dm_c,'file')
        createDesignMatrix(subject,nPreds,cond_locs_name{run},dir_dm_b,dir_dm_c);
    end
    
    dm_block{run} = dir_dm_b;       % Design matrix (single trial or single block wise)
    dm_cond{run} = dir_dm_c;        % Design matrix (condition wise)
    
    % POI: Patch-of-Interest (or ROI: Regresion-of-Interest)
    if POIfile_ind == 1
        poi_name = [dir_rois,'Auditory.nii']; % Auditory cortex
    elseif POIfile_ind == 2
        poi_name = [dir_rois, 'Motor.nii']; % Motor cortex
    else
        poi_name = [dir_rois,'Visual.nii']; % Visual cortex
    end
    
end

% MATLAB object definition %
input_file_info.dir_name = [dir_name,'sub-0',int2str(subject),'/'];
input_file_info.poi_name = poi_name;
input_file_info.func_name = func_name;
input_file_info.dm_name = dm_block;
input_file_info.dm2 = dm_cond;
input_file_info.CondClass = CondClass;
input_file_info.subject = subject;
input_file_info.cond_locs = cond_locs_name;
input_file_info.Hem = Hem;

pars(1) = nVols;
pars(2) = nPreds;
pars(3) = nTrials;
pars(4) = nPerRun;

input_file_info.pars = pars;
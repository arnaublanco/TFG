% Function that obtains the information from files.
%  INPUT:
%   · subject: Subject number
%   · Hem: 'LH' or 'RH'
%   · CondClass: Conditions to classify
%   · POIfile_ind: Index in POI file (1 -> Auditory, 2 -> Motor or 3 -> EVC) 
%  OUTPUT:
%   · input_file_info: MATLAB object containing the data info.
    
function [input_file_info] = getFileInfo(subject, Hem, CondClass, POIfile_ind, dataType)

if dataType == 1
    dir_name = '/Users/blancoarnau/Google Drive/TFG Arnau Blanco/data/output/'; % Path to the files
    dir_rois = 'ROIs/'; % Path to the ROIs
else
    dir_name = '/Users/blancoarnau/Documents/ds002715/derivatives/fmriprep/'; % Path to the files
    dir_rois = '/Users/blancoarnau/Documents/ds002715/derivatives/freesurfer/'; % Path to the ROIs/POIs
    dir_events = '/Users/blancoarnau/Documents/ds002715/'; % Path to events.tsv
end
    
% Fewer volumes are selected because there is a glitch at the beginning of the signal.
if subject == 1
   nVols = 113;
else
    nVols = 218;
end
nPreds = 18+1;  % Number of stimulation blocks plus 1 for baseline 
nTrials = 18;   % Number of stimulation blocks (without baseline) per run
nPerRun = 6;    % Number of blocks per condition per run
nRuns = 4;      % Number of runs

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

for run = 1:nRuns
    
    if dataType == 1
        % Functional image file
        func_name{run} = [dir_name,'sub-0',int2str(subject),'/sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-',int2str(run),'_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'];
    
        % File that maps single trials or blocks onto experimental conditions.
        cond_locs_name{run} = [dir_name,'sub-0',int2str(subject),'/sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-0',int2str(run),'_events.tsv'];
    else
        % Functional image file
        file = ['sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-',int2str(run),'_space-fsnative_hemi-',h,'_bold.func.gii'];
        func_name{run} = [dir_name,'sub-0',int2str(subject),'/ses-mri/func/',file];
        
        % Anatomical image file
        anat_name = [dir_name,'sub-0',int2str(subject),'/ses-mri/anat/','sub-0',int2str(subject),'_ses-mri_run-1_hemi-',h,'_inflated.surf.gii'];
        
        % File that maps single trials or blocks onto experimental conditions.
        cond_locs_name{run} = [dir_events,'sub-0',int2str(subject),'/ses-mri/func/sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-0',int2str(run),'_events.tsv'];
    end
    % Directories for the design matrices (single block-wise and condition-wise)
    dir_dm_b = [dir_name,'sub-0',int2str(subject),'/sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-0',int2str(run),'_block.xls'];
    dir_dm_c = [dir_name,'sub-0',int2str(subject),'/sub-0',int2str(subject),'_ses-mri_task-AVScenes',w,'_run-0',int2str(run),'_cond.xls'];
    
    % If design matrices files already exist, remove them.
    if exist(dir_dm_b,'file')
        delete(dir_dm_b);
    end
    
    if exist(dir_dm_c,'file')
        delete(dir_dm_c);
    end
    
    % Create design matrix and save them as .xls files
    createDesignMatrix(subject,nPreds,cond_locs_name{run},dir_dm_b,dir_dm_c);
    
    dm_block{run} = dir_dm_b;       % Design matrix (single trial or single block wise)
    dm_cond{run} = dir_dm_c;        % Design matrix (condition wise)
    
    % POI: Patch-of-Interest (or ROI: Region-of-Interest)
    if dataType == 1
        if POIfile_ind == 1
            poi_name = [dir_rois,'Auditory.nii']; % Auditory cortex
        elseif POIfile_ind == 2
            poi_name = [dir_rois, 'Motor.nii']; % Motor cortex
        else
            poi_name = [dir_rois,'Visual.nii']; % Visual cortex
        end
    else
        if POIfile_ind == 1
            poi_name = [dir_rois,'sub-0',int2str(subject),'/label/',lower(Hem),'.aparc.a2009s.annot']; % Auditory (tranverse temporal Gyrus -> Destrieux atlas)
        elseif POIfile_ind == 2
            %poi_name = [dir_rois,'sub-0',int2str(subject),'/label/',lower(Hem),'.aparc.a2009s.annot'];
            poi_name = [dir_rois,'sub-0',int2str(subject),'/label/',lower(Hem),'.BA4a_exvivo.label']; % Motor (primary motor area -> Brodmann area 4)
        else
            poi_name = [{[dir_rois,'sub-0',int2str(subject),'/label/',lower(Hem),'.V1_exvivo.label']},{[dir_rois,'sub-0',int2str(subject),'/label/',lower(Hem),'.V2_exvivo.label']}]; % Visual (V1 and V2)
        end
    end
    
end

% MATLAB object definition
input_file_info.dir_name = [dir_name,'sub-0',int2str(subject),'/'];
input_file_info.poi_name = poi_name;
input_file_info.func_name = func_name;
input_file_info.dm_name = dm_block;
input_file_info.dm2 = dm_cond;
input_file_info.CondClass = CondClass;
input_file_info.subject = subject;
input_file_info.cond_locs = cond_locs_name;
input_file_info.Hem = Hem;
if dataType == 2
    input_file_info.anat_name = anat_name;
end

pars(1) = nVols;
pars(2) = nPreds;
pars(3) = nTrials;
pars(4) = nPerRun;

input_file_info.pars = pars;
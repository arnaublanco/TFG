% Function that loads the specified Mesh Timecourse with a specified mask
% and hemisphere.
%  INPUT:
%   · subject: Subject number
%   · Patch_ind: Patch within the POIs that will be used (1, 2, 3, ... or 16) 
%   · Hem: 'LH' or 'RH'
%   · CondClass: Conditions to classify
%   · POIfile_ind: Index in POI file (1 -> Auditory, 2 -> Motor or 3 -> EVC) 
%   · visualize: Visualization of files when equal to 1.
%  OUTPUT:
%   · outD: MATLAB object containing the data.

function [outD] = readDataMtcPoi(subject, Patch_ind, Hem, CondClass, POIfile_ind, visualize)

addpath(genpath('/Users/blancoarnau/Documents/GitHub/gifti')) % Add GifTI library to import files
addpath(genpath('/Applications/freesurfer/matlab'));

fileNames = getFileInfo(subject, Hem, CondClass, POIfile_ind);  % Import file for a given hemisphere, subject and task

zBetas = 0;           % Z-score betas AFTER GLM prior to classifier analysis

zTimeSeries = 1;      % If zTimeSeries = 1, then the data is z-scored

% Sort file names from input
dirName = fileNames.dir_name;
poiName = fileNames.poi_name;
funcName = fileNames.func_name;
anatName = fileNames.anat_name;
dmName = fileNames.dm_name;
condLocs = fileNames.cond_locs;
p = fileNames.pars;  % nClass and nVols
subject = fileNames.subject;

% Parameters
nRuns = length(funcName);
nVols = p(1);
nPreds = p(2);
nTrials = p(3);
nPerRun = p(4);
CondClass = fileNames.CondClass;
nConditions = length(CondClass);

% Load POI file
if(POIfile_ind ~= 3)
    [l] = read_label('',poiName); % Read .label file
    locsV = l(:,1) + 1; % Save vertices in locsV (plus one because in MATLAB indices start at 1)
else
    if(exist('Patch_ind', 'var'))
        [l] = read_label('',poiName{Patch_ind}); % Read .label file
        locsV = l(:,1)+1; % Save vertices in locsV - this will be modified later - i want the mask to have 1s', 2s', 3s' etc
        %locsV = niftiread(find(poiName == Patch_ind));
    end
end

nVert = length(locsV);
funcData = zeros(nVols,nVert,nRuns);
DM = zeros(nVols,nPreds,nRuns);
betas = zeros(nPreds - 1,nVert,nRuns);  % Minus one because of mean confound
tvals = zeros(size(betas));

% In the following, a GLM is run for each run with the design matrix of each run.
% This basically does the "averaging" across time (functional data is a time series).
% The resulting beta weights for each voxel (derived for each block)
% are then used as input for the pattern classification, i.e. there will be
% a multivariate pattern of betas across voxels (space) but not across time.

% Main loop (across runs)
for r = 1:nRuns

    % Load Functional data
    main = gifti(convertStringsToChars(funcName{r})); % NOTE: Gifti does not work with strings -> conversion to characters
    
    main_masked = main.cdata(locsV,end-nVols+1:end); % Crop functional data with nVols and get vertices within POI    
    funcData(:,:,r) = main_masked';  % Save functional data

%     if(~isempty(find(funcData(:,:,r)==0)))
%         error('Zeros in Functional data: requires further thought!');
%     end

    % Visualize surfaces if set to 1
    close all;
    if visualize == 1
        anat = gifti(convertStringsToChars(anatName)); % Load anatomical file
        bold_masked = main; % Copy functional file to ´bold_masked´
        bold_masked.cdata(:,:) = 0; % Set all the vertices to 0
        bold_masked.cdata(locsV,1) = 1; % Set the vertices in the ROI to 1
        
        figure;
        plot(anat,main); colorbar; % Show functional file at t = 0
        title("Subject " + string(subject) + " - Run " + string(r) + " - " + Hem);
        
        figure;
        plot(anat,bold_masked); % Show mask
        title('Mask ' + string(POIfile_ind) + " - " + Hem);
        
        pause(5); % Pause for 5 seconds
    end
    
    % Load in design matrix file
    dmTmp = readDM(dmName{r},nPreds);
    DM(:,:,r) = dmTmp(end-nVols+1:end,:);

    if(size(DM,1)~=size(funcData,1))
        error('Design Matrix and functional data volume number does not match');
    end

    % PERFORM GLM COMPUTATION - SINGLE TRIAL / BLOCK COMPUTATION
    
    % GLM performed single-block, not deconvolved (can't be)
    % The timeseries are z-scored before running GLM in order to obtain comparable values to output
    [out,out2] = computeGLM(funcData(:,:,r), DM(:,:,r), zTimeSeries,zBetas);

    betas(:,:,r) = out(1:nPreds-1,:);  % Remove last beta, mean confound
    tvals(:,:,r) = out2(1:nPreds-1,:);  % t-values
    
end
% End loop across runs
    
% GET THE DESIGN SEQUENCE AND PARSE BETAS CONDITION-WISE

% Load sequence (condition data)
nTrials = size(betas,1);
%nPerRun = 6;
betasC = zeros(nPerRun,nVert,nConditions,nRuns);  % Betas split by condition
tvalsC = zeros(nPerRun,nVert,nConditions,nRuns);

for r = 1:nRuns
    
    %seq = condLocs(:,r);
    seq = readtsv(condLocs{r});
    seq = seq{3}(~strcmp(seq{3},{'baseline'}));
    %seq=seq';  %% flip for my data (in rows)
    %seq=seq(:);  %%concatenate into one trial vector of conditions codes

    for i = 1:nConditions
        f = [];
        for j = 1:size(seq,1)
            % Find where each condition is present (across trials)
            if strcmp(seq{j},CondClass{i}) 
                f = cat(2,f,j);
            end
        end
        % Separate out each condition
        if(length(f) == nPerRun)
            betasC(:,:,i,r) = betas(f,:,r);
            tvalsC(:,:,i,r) = tvals(f,:,r);
        else
            error('Not correct number per condition');
        end
    end
end

% Put data in structs
p(1) = nTrials;
p(2) = nConditions;
p(3) = nPerRun;
p(4) = nVert;  
p(5) = nRuns;
p(6) = nVols;

s = cell(1,4);
s{1,1} = subject;
s{1,2} = dirName;         % Where to save any output files
s{1,3} = p;               % Useful parameters
s{1,4} = locsV;           % All vertices in POI
s{1,5} = fileNames;   

outD = [];
outD.betasC = betasC;     % Betas by condition
outD.betas = betas;       % Betas by runs
outD.tvalsC = tvalsC;     % t-vals by condition
outD.data = funcData;     % Raw vertice timecourses for POI
outD.DM = DM;             % Design matrices
outD.S = s;               % Useful parameters

%%%% FINALLY COMPUTE UNIVARIATE ANOVAS  
% ** be careful ** not to select voxels on this basis for entry to classifier -
% is biased if selected across all runs, needs to be only training runs!!!
%anova = compute_ANOVAS(outD);  % compute univariate discrimination across ALL RUNS
%outD.anova = anova;
% this univariate analysis is run just to see whether there are univariate
% effects in the data. This is independent from the pattern classification.






function [outD] = readDataMtcPoi(subject, Patch_ind, Hem, CondClass, POIfile_ind, dataType)

% INPUT:
    % subject: Subject number
    % Patch_ind: Patch within the POIs that will be used (1, 2, 3, ... or 16) 
    % Hem: 'LH' or 'RH'
    % CondClass: Conditions to classify
    % POIfile_ind: Index in POI file (1 -> Auditory, 2 -> Motor or 3 -> EVC)
    % dataType: Volume-based data (1) or surface-based data (2).
    
% OUTPUT:
    % outD: MATLAB object containing the data.

fileNames = getFileInfo(subject,Hem,CondClass,POIfile_ind,dataType);  % Import file for a given hemisphere, subject and task

zBetas = 0;           % Z-score betas AFTER GLM prior to classifier analysis

zTimeSeries = 1;      % If zTimeSeries = 1, then the data is z-scored

% Get filenames data from input
dirName = fileNames.dir_name;
poiName = fileNames.poi_name;
funcName = fileNames.func_name;
dmName = fileNames.dm_name;
condLocs = fileNames.cond_locs;
p = fileNames.pars;
subject = fileNames.subject;

% Parameters
nRuns = length(funcName);
nVols = p(1);
nPreds = p(2);
nPerRun = p(4);
CondClass = fileNames.CondClass;
nConditions = length(CondClass);

% Load in POI/ROI file
if dataType == 1
    locsV = niftiread(poiName);
    
    idx = find(locsV == Patch_ind);
    [x,y,z] = ind2sub(size(locsV),idx);

    if(strcmp(Hem,'LH'))
        pos = x < size(locsV,1)/2;
        x = x(pos);
        y = y(pos);
        z = z(pos);
    else
        pos = x > size(locsV,1)/2;
        x = x(pos);
        y = y(pos);
        z = z(pos);
    end

    nVox = sum(pos,'all');
else
    if(POIfile_ind == 1)
        [~, label, ~] = read_annotation(poiName); % Read .annot file
        locsV = find(label == 14433340); % Tranverse temporal gyrus has label 14433340
    elseif(POIfile_ind == 2)
        [l] = read_label('',poiName); % Read .label file
        locsV = l(:,1) + 1; % Save vertices in locsV (plus one because in MATLAB indices start at 1)
    elseif(POIfile_ind == 3)
        if(exist('Patch_ind', 'var'))
            [l] = read_label('',poiName{Patch_ind}); % Read .label file
            locsV = l(:,1) + 1; % Save vertices in locsV (plus one because in MATLAB indices start at 1)
        end
    end
    nVox = length(locsV);
end

funcData = zeros(nVols,nVox,nRuns); % Functional data
DM = zeros(nVols,nPreds,nRuns); % Design matrix
betas = zeros(nPreds-1,nVox,nRuns); % Betas across blocks
tvals = zeros(size(betas)); % T-values

% In the following, a GLM is performed for each run with their corresponding design matrix.
% This basically does an "average" across time, as functional data is a time series.
% The resulting beta weights are for each voxel and are used as input for
% the pattern classification, i.e. betas across voxels (space) but not across time.

% Main loop - Data import and GLM
fprintf('\n- Loading fMRI scans...');
for r = 1:nRuns

    if dataType == 1
        % Load in functional data
        main = niftiread(funcName{r});
        main = main(:,:,:,end-nVols+1:end); % Crop to nVols

        % Select voxels in the ROI
        tmp = zeros(nVox,nVols);
        for i = 1:size(x,1)
            tmp(i,:) = reshape(main(x(i),y(i),z(i),:),[1 nVols]);
        end

        % Tranpose to fit in funcData
        funcData(:,:,r) = tmp';
    else
        % Load functional data
        main = gifti(convertStringsToChars(funcName{r})); % NOTE: Gifti does not work with strings -> conversion to characters
    
        tmp = main.cdata(locsV,end-nVols+1:end); % Crop functional data with nVols and get vertices within POI
    
        funcData(:,:,r) = tmp';  % Save functional data
    end
    
    % Load in Design Matrix
    dmTmp = readmatrix(dmName{r});
    DM(:,:,r) = dmTmp(end-nVols+1:end,:);

    if(size(DM,1)~=size(funcData,1))
        error('Design Matrix and functional data volume number does not match');
    end
    
    % GLM computation (single block/trial computation)
    fprintf('\n- Computing betas for run '+string(r)+'...');
    [out,out2] = computeGLM(funcData(:,:,r), DM(:,:,r), zTimeSeries,zBetas);

    betas(:,:,r) = out(2:end,:);   % Remove first beta (dummy variables' weight)
    tvals(:,:,r) = out2(2:end,:);  % Remove first t-value
    
end

% Load in design sequence (condition data)
nTrials = size(betas,1);
betasC = zeros(nPerRun,nVox,nConditions,nRuns);  % Betas split by condition
tvalsC = zeros(nPerRun,nVox,nConditions,nRuns);  % Betas split by block/trial

for r = 1:nRuns
    
    seq = readtsv(condLocs{r}); % Read design sequence file
    seq = seq{3}(~strcmp(seq{3},{'baseline'}));

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

% Format the data for outputting
p(1) = nTrials;
p(2) = nConditions;
p(3) = nPerRun;
p(4) = nVox;  
p(5) = nRuns;
p(6) = nVols;

s = cell(1,4);
s{1,1} = subject;
s{1,2} = dirName;
s{1,3} = p;
s{1,4} = locsV;
s{1,5} = fileNames;   

outD = [];
outD.betasC = betasC;
outD.betas = betas;
outD.tvalsC = tvalsC;
outD.data = funcData;
outD.DM = DM;
outD.S = s;

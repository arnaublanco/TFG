function [outD] = readDataMtcPoi(subject, Patch_ind, Hem, CondClass, POIfile_ind)

% INPUT:
    % subject: Subject number
    % Patch_ind: Patch within the POIs that will be used (1, 2, 3, ... or 16) 
    % Hem: 'LH' or 'RH'
    % CondClass: Conditions to classify
    % POIfile_ind: Index in POI file (1 -> Auditory, 2 -> Motor or 3 -> EVC)
    
% OUTPUT:
    % outD: MATLAB object containing the data.

fileNames = getFileInfo(subject,Hem,CondClass,POIfile_ind);  % Import file for a given hemisphere, subject and task

zBetas = 0;           % Z-score betas AFTER GLM prior to classifier analysis

zTimeSeries = 1;      % If zTimeSeries = 1, then the data is z-scored

%%% Sort file names from input
dirName = fileNames.dir_name;
poiName = fileNames.poi_name;
funcName = fileNames.func_name;
dmName = fileNames.dm_name;
condLocs = fileNames.cond_locs;
p = fileNames.pars;  %% nClass and nVols
subject = fileNames.subject;

%%% Parameters
nRuns = length(funcName);
nVols = p(1);
nPreds = p(2);
nTrials = p(3);
nPerRun = p(4);
CondClass = fileNames.CondClass;
nConditions = length(CondClass);

% Load POI file
locsV = [];
if(POIfile_ind ~= 3)
    locsV = niftiread(poiName);       
    %they have used ´unique´ here, which returns the same array without repetitions
else
    if(exist('Patch_ind', 'var'))
        locsV = niftiread(poiName); % this will be modified later - i want the mask to have 1s', 2s', 3s' etc
        %locsV = niftiread(find(poiName == Patch_ind));
    end
end

idx = find(locsV == 1);
[x,y,z] = ind2sub(size(locsV),idx);

if(Hem == 'LH')
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
funcData = zeros(nVols,nVox,nRuns);
DM = zeros(nVols,nPreds,nRuns);
betas = zeros(nPreds-1,nVox,nRuns);  %% minus one for mean confound
tvals = zeros(size(betas));

% In the following, a GLM is run for each run with the design matrix of each run.
% This basically does the "averaging" across time (functional data is a time series).
% The resulting beta weights for each voxel (derived for each block)
% are then used as input for the pattern classification, i.e. there will be
% a multivariate pattern of betas across voxels (space) but not across time.

% Main loop
fprintf('\n- Loading fMRI scans...');
for r = 1:nRuns

    % Load in Functional data
    main = niftiread(funcName{r});
    main = main(:,:,:,end-nVols+1:end);
    
    tmp = zeros(nVox,nVols);
    for i = 1:size(x,1)
        tmp(i,:) = reshape(main(x(i),y(i),z(i),:),[1 nVols]);
    end
    
    funcData(:,:,r) = tmp';
    
    % Load in DESIGN MATRIX file
    dmTmp = readmatrix(dmName{r});
    DM(:,:,r) = dmTmp(end-nVols+1:end,:);

    if(size(DM,1)~=size(funcData,1))
        error('Design Matrix and functional data volume number does not match');
    end

    % PERFORM GLM COMPUTATION - SINGLE TRIAL / BLOCK COMPUTATION
    
    % GLM performed single trial / block, not deconvolved (can't be)
    % importantly I am z-scoring the timeseries here before running GLM
    % in order to obtain comparable values to output
    fprintf('\n- Computing betas for run '+string(r)+'...');
    [out,out2] = computeGLM(funcData(:,:,r), DM(:,:,r), zTimeSeries,zBetas);

    betas(:,:,r) = out(1:nPreds-1,:);  %% remove last beta, mean confound
    tvals(:,:,r) = out2(1:nPreds-1,:);  %% t values
    %the dimensions here (e.g. 1:nPreds-1) need to match whatever was defined above,
    %otherwise dimension mismatch
    
end   %% end loop across runs
    
% GET THE DESIGN SEQUENCE AND PARSE BETAS CONDITION-WISE

% load in sequence (condition data)
nTrials = size(betas,1);
%nPerRun=6;
betasC = zeros(nPerRun,nVox,nConditions,nRuns);  %% betas split by condition
tvalsC = zeros(nPerRun,nVox,nConditions,nRuns);

for r = 1:nRuns
    
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


%PUT IN NICE STRUCTURES for OUTPUTTING
p(1) = nTrials;
p(2) = nConditions;
p(3) = nPerRun;
p(4) = nVox;  
p(5) = nRuns;
p(6) = nVols;

s = cell(1,4);
s{1,1} = subject;
s{1,2} = dirName;         %% where to save any output files
s{1,3} = p;               %% useful parameters
s{1,4} = locsV;           %% all verts in POI
s{1,5} = fileNames;   


outD = [];
outD.betasC = betasC;     %% betas by condition
outD.betas = betas;       %% betas by runs
outD.tvalsC = tvalsC;     %% t vals by condition
outD.data = funcData;      %% raw vertice timecourses for POI
outD.DM = DM;             %% design matrices
outD.S = s;               %% useful parameters



%%%% FINALLY COMPUTE UNIVARIATE ANOVAS  
% ** be careful ** not to select voxels on this basis for entry to classifier -
% is biased if selected across all runs, needs to be only training runs!!!
%anova=compute_ANOVAS(outD);  %% compute univariate discrimination across ALL RUNS
%outD.anova=anova;
%this univariate analysis is run just to see whether there are univariate
%effects in the data. This is independent from the pattern classification.
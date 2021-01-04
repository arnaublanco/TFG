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

if(Hem == 'LH')
    locsV = locsV(1:floor(size(locsV,1)/2),:,:);
else
    locsV = locsV(ceil(size(locsV,1)/2):end,:,:);
end

nVox = sum(locsV,'all');
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
for r = 1:nRuns

    % Load Functional data
    main = niftiread(funcName{r});
    
    idx = find(repmat(locsV,[1 1 1 nVols])); % Find indices where there's a 1.
    main_masked = main(:,:,:,end-nVols+1:end); % Crop functional data with nVols
    
    funcData(:,:,r) = reshape(main_masked(idx),[nVox nVols])';  %% get timecourses from POI file

    if(~isempty(find(funcData(:,:,r)==0)))
        error('Zeros in Functional data: requires further thought!');
    end
    
    % Load in DESIGN MATRIX file
    dmTmp = readDM(dmName{r},nPreds);
    DM(:,:,r) = dmTmp(end-nVols+1:end,:);

    if(size(DM,1)~=size(funcData,1))
        error('Design Matrix and functional data volume number does not match');
    end

    % PERFORM GLM COMPUTATION - SINGLE TRIAL / BLOCK COMPUTATION
    
    % GLM performed single trial / block, not deconvolved (can't be)
    % importantly I am z-scoring the timeseries here before running GLM
    % in order to obtain comparable values to output
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
    
    %seq = condLocs(:,r);
    seq = readtsv(condLocs{r});
    seq = seq{3}(~strcmp(seq{3},{'baseline'}));
    %seq=seq';  %% flip for my data (in rows)
    %seq=seq(:);  %%concatenate into one trial vector of conditions codes

    for i = 1:nConditions
        f = [];
        for j = 1:size(seq,1)
            if strcmp(seq{j},CondClass{i}) 
                f = [f j];
            end %% find where each cond is present (across trials)
        end
        if(length(f) == nPerRun)
            betasC(:,:,i,r) = betas(f,:,r);   %% separate out each condition
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






function run_single_sub_CollapseHem_main_analysis_parallel(subject,Patch_ind,CondClass,POIfile_ind,visualize)

%%% Script for permutation analysis
% subject = 'CGY09';
% CondClass = [1 2 3];
% Patch_ind = 1;
% POIfile_ind = 2;

tic % Start a stopwatch timer.

%%% Script needs to correspond to the numbers (condition definitions) in the trial sequence
% subject='DHE05';
% Patch_ind=1;  %% poi index (eg. 1 for V1, 2 for V2v etc. depending on how
% many pois and in which order they are defined)
% CondClass=[1 2 3];   %% CRUCIAL- which conditions to classify

% CHANGE file name under which the results are saved - see below

permGP = 1; % 0: Standard analysis; 1: Permutation analysis

% Prepare data - get timeseries, do GLM etc
[outDL] = readDataMtcPoi(subject, Patch_ind, 'LH', CondClass, POIfile_ind, visualize); % Left hemisphere
[outDR] = readDataMtcPoi(subject, Patch_ind, 'RH', CondClass, POIfile_ind, visualize); % Right hemisphere

% Concatenate across vertice dimension - hemispheres
outD = [];               % Define empty variable outD
f1 = size(outDL.betasC); % Size of matrix of betas from left hemisphere
f2 = size(outDR.betasC); % Size of matrix of betas from right hemisphere

% Check that dimensions match
if(f1(1) == f2(1) && f1(3) == f2(3) && f1(4) == f2(4))
    betasC = cat(2,outDL.betasC,outDR.betasC);  % Concatenate condition-wise betas matrices (left and right hemispheres)
    betas = cat(2,outDL.betas,outDR.betas);     % Concatenate block-wise betas matrices (left and right hemisphere)
    S{3} = outDL.S{3};                          % Permutation labels
    S{3}(4) = f1(2)+f2(2);                      % Sum of dimensions of both hemispheres 
    f = [];                                     % Initialize f
    f.CondClass = CondClass;                    % Add conditions
    S{5} = f;                                   % Save struct with conditions in S
    outD.betasC = betasC;                       % Save betas (condition-wise) in outD
    outD.betas = betas;                         % Save betas (block-wise) in outD
    outD.S = S;                                 % Save S in outD
    
else
    error('Hemispheric Data mismatch');
end

% Load the randomisation vector for the permutation analysis (same for
% every subject)
load randVecPerm;
inputRandVec = f;

% Perform classification for Permutation Analysis
Spc = zeros(1000,1);
Apc = zeros(1000,1);
betasC = outD.betasC;
p = outD.S{3};
CondClass = outD.S{5}.CondClass;

toc % Elapsed time for GLM, load the data, etc

tic % Start a stopwatch timer

% Parse for cross-validation cycles
[train_set, test_set, anovas] = parse_runs_surf(betasC);

toc % Elapsed time for parsing

tic % Start a stopwatch timer

% ´parfor´: Executes for loop in parallel with other processes %
nPermutations = 1000;
parfor (i = 1:nPermutations)
    % randVec = inputRandVec(:,i);
    [svmOutPerm] = singleSVM_PermP(train_set, test_set, p, CondClass, permGP,inputRandVec(:,i));
    Spc(i) = mean(svmOutPerm.pc);
    Apc(i) = mean(svmOutPerm.av);
end

toc  % Elapsed for time for 1000 permutations in parallel

% Perform classification for real data %
permGP = 0;
[svmOutObs] = singleSVMP(betasC, p, CondClass, permGP);
Obs_Spc = mean(svmOutObs.pc);
Obs_Apc = mean(svmOutObs.av);

pPerm_Spc = length(find(Spc >= Obs_Spc)) ./ 1000;
pPerm_Apc = length(find(Apc >= Obs_Apc)) ./ 1000;

% t test across runs
%nConditions=length(outD.S{5}.CondClass);  %% conditions being classified
%[hST,pST,ci,statsST]=ttest(svmOut.pc,1/nConditions,.05,'right'); % single trial
%[hAV,pAV,ci,statsAV]=ttest(svmOut.av,1/nConditions,.05,'right'); % average

%save FWStestCGY09_V5.mat
% to save the output
outname = sprintf('%s_MainAnalysis_CollapseHem_Patch%d_POI%d.mat', subject, Patch_ind,POIfile_ind);
save(outname, 'subject','Patch_ind','permGP','outD','svmOutObs','pPerm_Spc','pPerm_Apc','Spc','Apc');

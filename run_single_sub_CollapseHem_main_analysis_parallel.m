% Function that computes classification for a given subject and ROI.
%  INPUT:
%   · subject: Subject number.
%   · Patch_ind: Patch within ROI (only used for EVC).
%   · CondClass: Experimental conditions.
%   · POIfile_ind: POI/ROI (Visual, Motor or EVC).
%   · dataType: Volume-based (1) or surface-based (2).
%   · classifier: SVM (1) or SVM-RFE (2).
%   · cv: Cross-validation per run (1) or per block (2).

function run_single_sub_CollapseHem_main_analysis_parallel(subject, Patch_ind, CondClass, POIfile_ind, dataType, classifier, cv)

% tmp = subject;
% if POIfile_ind == 1
%     load('Results/SVM/surface/across_blocks/sub-01_MainAnalysis_CollapseHem_Patch1_POI1.mat')
% elseif POIfile_ind == 2
%     load('Results/SVM/surface/across_blocks/sub-01_MainAnalysis_CollapseHem_Patch1_POI2.mat')
% else
%    if Patch_ind == 1
%        load('Results/SVM/surface/across_blocks/sub-01_MainAnalysis_CollapseHem_Patch1_POI3.mat')
%    else
%        load('Results/SVM/surface/across_blocks/sub-01_MainAnalysis_CollapseHem_Patch2_POI3.mat')
%    end
% end
% subject = tmp;

if POIfile_ind == 1
    roi = 'auditory cortex';
elseif POIfile_ind == 2
    roi = 'motor cortex';
else
    if Patch_ind == 1
        roi = 'EVC - V1';
    else
        roi = 'EVC - V2';
    end
end

if classifier == 1
    folder = 'SVM';
else
    folder = 'SVM_RFE';
    cv = 2;
end

fprintf('-------------------------------------\n');
fprintf('Computing: sub-0'+string(subject)+' ('+string(roi)+')');
tic % Start a stopwatch timer.

% Prepare data - get timeseries, do GLM etc
fprintf('\n\nLeft hemisphere:');
[outDL] = readDataMtcPoi(subject, Patch_ind, 'LH', CondClass, POIfile_ind, dataType); % Left hemisphere
fprintf('\n\nRight hemisphere:');
[outDR] = readDataMtcPoi(subject, Patch_ind, 'RH', CondClass, POIfile_ind, dataType); % Right hemisphere

% Concatenate across vertice dimension - hemispheres
outD = [];               % Define empty variable outD
f1 = size(outDL.betasC); % Size of matrix of betas from left hemisphere
f2 = size(outDR.betasC); % Size of matrix of betas from right hemisphere

% Check that dimensions match
fprintf('\n\nJoining betas from LH and RH... ');
pause(1);
if(f1(1) == f2(1) && f1(3) == f2(3) && f1(4) == f2(4))
    betasC = cat(2,outDL.betasC,outDR.betasC);    % Concatenate block-wise betas matrices (left and right hemisphere)
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
if cv == 1
    load randVecPerm;
else
    load randVecperm_block;
end
inputRandVec = f;

% Perform classification for Permutation Analysis
nPermutations = 1000;
Spc = zeros(nPermutations,1);
Apc = zeros(nPermutations,1);
betasC = outD.betasC;
p = outD.S{3};
CondClass = outD.S{5}.CondClass;
nCond = 1:size(CondClass,2);

toc % Elapsed time for GLM, load the data, etc

tic % Start a stopwatch timer

% Perform classification for real data
fprintf(['\nComputing ',folder,' without permutation of labels...']);
pause(1);
permGP = 0; % 0: Standard analysis; 1: Permutation analysis
if classifier == 1
    if cv == 1
        [OutObs] = singleSVMP(betasC, nCond, permGP);
    else
        [OutObs] = singleSVMP_block(betasC, nCond, permGP);
    end
else
    [OutObs] = singleSVMP_RFE(betasC, nCond);
end

Obs_Spc = mean(OutObs.pc); % Percentage correct single-block
Obs_Apc = mean(OutObs.av); % Percentage correct average

toc % Elapsed time for SVM without permutation of labels
pause(1);
tic % Start a stopwatch timer

% Executes for loop in parallel with other processes
fprintf(['\nComputing ',folder,' with permuted labels...']);
pause(1);
permGP = 1;
% Parse for cross-validation cycles
if cv == 1
    [train_set, test_set] = parse_runs_surf(betasC);
else
    [train_set, test_set] = parse_runs_surf_blocks(betasC);
end

if classifier ~= 3
    parfor (i = 1:nPermutations)
        if cv == 1
            [OutPerm] = singleSVM_PermP(train_set, test_set, p, 1:3, permGP, inputRandVec(:,i));
        else
            [OutPerm] = singleSVMP_PermP_block(train_set, test_set, p, 1:3, permGP, inputRandVec(:,i));
        end
        Spc(i) = mean(OutPerm.pc);
        Apc(i) = mean(OutPerm.av);
    end
else
    parfor (i = 1:nPermutations)
        [OutPerm] = singleKNN_Perm(train_set, test_set, p, 1:3, permGP, inputRandVec(:,i));
        Spc(i) = mean(OutPerm.pc);
        Apc(i) = mean(OutPerm.av);
    end
end
fprintf(['\n',int2str(nPermutations),' realizations with permuted labels done.\n']);

% toc % Elapsed for time for 1000 permutations in parallel


pPerm_Spc = length(find(Spc >= Obs_Spc)) ./ nPermutations;
pPerm_Apc = length(find(Apc >= Obs_Apc)) ./ nPermutations;

% t-test across runs
nConditions = length(outD.S{5}.CondClass);  % Conditions being classified
[hST,pST,ci,statsST] = ttest(OutObs.pc,1/nConditions,.05,'right'); % Single trial
[hAV,pAV,ci,statsAV] = ttest(OutObs.av,1/nConditions,.05,'right'); % Average

% Save FWStestCGY09_V5.mat to save the output
pause(1);
disp(['Saving results of ',folder,'...']);
outname = sprintf(['Results/',folder,'/sub-0%i_MainAnalysis_CollapseHem_Patch%d_POI%d.mat'], subject, Patch_ind, POIfile_ind);
save(outname, 'subject','Patch_ind','permGP','outD','OutObs','pPerm_Spc','pPerm_Apc','Spc','Apc');
end
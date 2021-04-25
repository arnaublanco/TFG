% Treball de Fi de Grau: Decoding Natural Sounds in "Early Visual" Cortex of Congenitally Blind Individuals.

% The code was originally created by Fraser W. Smith and was adapted by Petra Vetter and Lukasz Bola for the study of the paper.
% It has also been modified for my study by me, Arnau Blanco Borrego.

% Requirements: LibSVM and NeuroElf.
addpath('/Users/blancoarnau/Documents/MATLAB/libsvm/matlab/')
addpath('/Users/blancoarnau/Documents/GitHub/neuroelf-matlab/');

% Run run_single_sub_CollapseHem_main_analysis_parallel.m:
%   1. Preparation of input files: getFileInfo.m
%   2. Create Design Matrix: createDesignMatrix.m
%   3. Data import: readDataMtcPoi.m
%   4. Compute betas with GLM: computeGLM.m
%   5. Classification for every subject with and without permutation of labels.
run run_all_sub_main_analysis.m

% Run single_sub_meanClassAccu.m:
% Saves averaged accuracy for each subject and ROI across cross-classification folds. 
run single_sub_meanClassAccu.m

% Run groupAnalysis_Perm.m:
%   1. Compare actual classification accuracies with random permutations.
%   2. Obtain p-values.
run groupAnalysis_Perm.m

% Run multipleCompCorr_SingleThreshold.m:
%   1. Correct general analyses using single threshold test.
%   2. Two separate corrections are calculated - one for EVC, auditory cortex
%   and the motor cortex; and one for V1, V2 and V3
%run multipleCompCorr_SingleThreshold.m

% Correct analysis by eccentricity using FDR (fdr_bh function) - see the script for the list of visual areas included 
%run multipleCompCorr_fdr.m

%%% Run additional pipeline to get unbiased confusion matrices
% Perform SVM classification (note that "permGP" in run_single_sub_CollapseHem_get_unbiased_conf_mats.m 
% must be set to 0, to get actual classification accuracies, not the accuracies of random
% permutations)
%run run_all_sub_get_unbiased_conf_mats.m

% get the confusion matrices
%run pooled_confusionMatrix_unbiased_visual.m

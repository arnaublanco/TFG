% Treball de Fi de Grau: Decoding Natural Sounds in Early "Visual" Cortex of Congenitally Blind Individuals.

% The code was originally created by Fraser W. Smith and was adapted by Petra Vetter and Lukasz Bola for the study of the paper.
% It has also been modified by Arnau Blanco Borrego.

% Requirements: LibSVM, NeuroElf, GifTi and FreeSurfer.
addpath(genpath('/Users/blancoarnau/Documents/MATLAB/libsvm/matlab/'));
addpath(genpath('/Users/blancoarnau/Documents/GitHub/neuroelf-matlab/'));
addpath(genpath('/Users/blancoarnau/Documents/GitHub/gifti'))
addpath(genpath('/Applications/freesurfer/matlab'));

% Options: Run classification on volume-based data (1) and/or on
% surface-based data (2)
dataType = 1;

% Options: Classifier -> 1: SVM and 2: SVM-RFE
classifier = 1;

% Options: In the case of SVM, run cross-validation per run (1) or per
% block (2).
cv = 2;

% Raise error if the parameters do not have correct values.
if (dataType ~= 1 && dataType ~= 2) || (classifier ~= 1 && classifier ~= 2 && classifier ~= 3) || (cv ~= 1 && cv ~= 2)
   error('The selected parameters are not correct.'); 
end

if dataType == 1   
    fprintf('** Volume-based data ** \n');
else
    fprintf('** Surface-based data ** \n');
end

% Run run_single_sub_CollapseHem_main_analysis_parallel.m:
%   1. Preparation of input files: getFileInfo.m
%   2. Create Design Matrix: createDesignMatrix.m
%   3. Data import: readDataMtcPoi.m
%   4. Compute betas with GLM: computeGLM.m
%   5. Classification for every subject with and without permutation of labels.
run_all_sub_main_analysis(dataType,classifier,cv);

% Run single_sub_meanClassAccu.m:
% Saves averaged accuracy for each subject and ROI across cross-classification runs. 
single_sub_meanClassAccu(classifier);

% Run groupAnalysis_Perm.m:
%   1. Compare actual classification accuracies with random permutations.
%   2. Obtain p-values.
groupAnalysis_Perm(classifier);

disp('Pipeline completed.');

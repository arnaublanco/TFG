% Function that computes SVM with a permutation analysis that uses always
% the same randomization vector.
%  INPUT:
%   · train_set: Train set.
%   · test_set: Test set.
%   · p: Data info.
%   · CondClass: Conditions to classify (stimuli).
%   · permGP: Permutation type (0: normal analysis; 1: randomization of labels).
%   · inputRandVec: Randomization vector.
%  OUTPUT:
%   · svmOut: MATLAB object containing the results of the SVM.

function [svmOut] = singleSVM_PermP(train_set, test_set, p, CondClass, permGP, inputRandVec)

nConditions = p(2);
nPerRun = p(3);
nRuns = p(5);

% Output variables
svm_ws = cell(1,nRuns);
svm_class = zeros(nPerRun*(nConditions),nRuns);
svma_class = zeros(nConditions,nRuns);
svm_cm = zeros(nConditions,nConditions,nRuns);
svm_pc = zeros(nRuns,1);
svm_av = zeros(nRuns,1);
svm_mod = cell(1,nRuns);
   
% Main loop - Cross-validation cycles in every run  

for r = 1:size(train_set,3)  

    % Define train and test for this cycle of cross-validation
    train = train_set(:,:,r);  
    test2 = test_set(:,:,r);

    % Define coding variable for training - assumes LOOCV (leave one run out)   
    gp = []; k = 1; l = nPerRun*(nRuns-1);
    for ii = 1:nConditions
        gp(k:l,1) = ii;
        k = k + nPerRun*(nRuns-1);
        l = l + nPerRun*(nRuns-1);
    end
    
    % Permutation of label vector - this will allow the creation of a randomization
    % distribution for good comparison purposes.
    if(permGP == 1)
        % Not good for group analysis, as here the randomization vector must
        % be the same across subjects
        f = inputRandVec;
        gp = gp(f); % Take the (constant) randomization vector from the input
    end
    
    % If the training set and the labeling set do not have the dimensions,
    % then show an error.
    if(size(train,1) ~= size(gp,1))
        error('Training and gp vector mismatch');
    end

    % Run it on the average for each condition in test run
    in = zeros(nConditions, size(test2,2));
    k = 1; l = nPerRun;
    for i = 1:nConditions
        in(i,:) = mean(test2(k:l,:));
        k = k + nPerRun; l = l + nPerRun;
    end

    % *** Single trial SVM prediction ***
    
    % Set train set to -1 to 1 scale
    [train, pars] = stretch_cols_ind(train, -1, 1); 
    
    % Train SVM model
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,unique(v)) = 0; % Remove NaNs from train set
    
    svm_model = svmtrain(gp, train, '-t 0 -c 1'); % -t 0 = linear SVM, -c 1 = cost value of 1
    svm_mod{r} = svm_model;
    
    % Get the weights from model
    svm_weights = svm_DefineWeights(svm_model);  
    svm_ws{r} = svm_weights;

    % Define coding variable for testing
    gp_test = []; k = 1; l = nPerRun;
    for ii = 1:nConditions
        gp_test(k:l,1) = ii;
        k = k + nPerRun;
        l = l + nPerRun;
    end
    
    % Test SVM model
    
    % Set train set to -1 to 1 scale
    [test2] = stretchWithGivenPars(test2, [-1 1], pars);
    test2(:,unique(v)) = 0; % Remove NaNs from test set
    
    [svm_class(:,r), ~, ~] = svmpredict(gp_test, test2, svm_model);
    
    % Compute performance on testing runs
    f = zeros(nConditions);  % Confusion matrix 

    k = 1; l = nPerRun;
    for i = 1:nConditions
        for j = 1:nConditions
            f(i,j) = length(find(svm_class(k:l,r) == j));
        end
        k = k + nPerRun; l = l + nPerRun;
    end

    svm_cm(:,:,r) = f;
    svm_pc(r) = trace(f) / (nPerRun*(nConditions));  

    % *** Test SVM Prediction of Condition Average in Test Run ***
    [in] = stretchWithGivenPars(in, [-1 1], pars);
    in(:,unique(v)) = []; % Remove NaNs from test set in
    [svma_class(:,r), ~, ~] = svmpredict(CondClass', in, svm_model);

    svm_av(r) = length(find(svma_class(:,r) == CondClass')) ./ (nConditions);

end
% End cross-validation loop

% Proper output format
data{1} = train_set;
data{2} = test_set;

svmOut = [];
svmOut.models = svm_mod;    % SVM model for each cross-validation fold
svmOut.ws = svm_ws;         % Weights for each binary classification pbm
svmOut.class = svm_class;   % Single-block classifications
svmOut.Aclass = svma_class; % Average classification
svmOut.cm = svm_cm;         % Confusion matrices for single block classifications
svmOut.pc = svm_pc;         % Percentage correct single-block
svmOut.av = svm_av;         % Percentage correct average
svmOut.data = data;         % Parsed data used in classifier

% t test on performance across runs, within subject
[h,p,ci,stats] = ttest(svmOut.pc,1/nConditions,.05,'right');
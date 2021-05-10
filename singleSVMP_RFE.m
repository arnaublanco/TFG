% Function that computes Support Vector Machine for the betas.
%  INPUT:
%   · betasC: Matrix containing the betas.
%   · CondClass: Stimuli (1: Forest, 2: People, 3: Traffic)
%  OUTPUT:
%   · svmOut: MATLAB object containing the results of the SVM.

function [svmOut] = singleSVMP_RFE(betasC, CondClass)

nConditions = size(betasC,3);
nPerRun = size(betasC,1);
nRuns = size(betasC,4);

fprintf('\nParsing train and test sets...\n');
% Parse for cross-validation cycles
[train_set, test_set, anovas] = parse_runs_surf_blocks(betasC);

% Output variables
nFolds = size(train_set,3);
svma_class = zeros(nConditions,nFolds);
svm_av = zeros(nFolds,1);
svm_mod = cell(1,nFolds);

% Define coding variable for training - assumes LOOCV (leaves one run out)
gp = []; k = 1; l = nPerRun*nRuns-1;
for ii = 1:nConditions
    gp(k:l,1) = ii;
    k = k + nPerRun*nRuns-1;
    l = l + nPerRun*nRuns-1;
end

% Define coding variable for testing
gp_test = []; k = 1; l = 1;
for ii = 1:nConditions
    gp_test(k:l,1) = ii;
    k = k + 1;
    l = l + 1;
end

acc_avg = []; % Initialize accuracies
ranks = []; % Initialize ranks
nMax = size(train_set,2); % Maximum number of ranks
if nMax > size(train_set,2)
    nMax = size(train_set,2);
end

% Recursive Feature Elimination: Cross-validation cycles in every run
for f = 1:size(train_set,3)

    % Define train set and test set for this cycle of cross-validation
    train = train_set(:,:,f);
    test2 = test_set(:,:,f);
    
    % If the training set and the labeling set do not have the dimensions,
    % then show an error.
    if(size(train,1) ~= size(gp,1))
        error('Training and gp vector mismatch');
    end

    % *** Single trial SVM prediction ***
    
    % Set train set to -1 to 1 scale
    [train, pars] = stretch_cols_ind(train, -1, 1); 
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,unique(v)) = 0; % Remove NaNs from train set
    
    % Set train set to -1 to 1 scale
    [test2] = stretchWithGivenPars(test2, [-1 1], pars);
    test2(:,unique(v)) = 0; % Remove NaNs from test set
    
    fprintf('\nTraining SVM model for fold ' + string(f) + ' with RFE...\n');
    pause(1);
    
    [ftRank,~] = ftSel_SVMRFECBR(train,gp); % Compute RFE
    acc = [];
    % Loop to compute accuracies with respect to the value of k
    for k = 1:nMax
        idx = find(ftRank <= k);
        train_k = train(:,idx);
        svm_model_k = svmtrain(gp, train_k, '-t 0 -c 1');
        [svm_class_k, ~, ~] = svmpredict(CondClass', test2, svm_model_k);
        acc_k = length(find(svm_class_k == CondClass')) ./ (nConditions);
        acc = cat(1,acc,acc_k(1));
    end
    acc_avg = cat(2,acc_avg,acc);
    ranks = cat(2,ranks,ftRank');
    
end

acc_avg = mean(acc_avg,2);
maxAcc = max(acc_avg);
numK = find(round(acc_avg,3) == round(maxAcc,3));
numK = numK(end);

% Main loop - Cross-validation cycles in every run
for f = 1:size(train_set,3)  

    % Define train set and test set for this cycle of cross-validation
    train = train_set(:,:,f);  
    test2 = test_set(:,:,f);
    
    % If the training set and the labeling set do not have the dimensions,
    % then show an error.
    if(size(train,1) ~= size(gp,1))
        error('Training and gp vector mismatch');
    end

    % *** Single trial SVM prediction ***
    
    % Set train set to -1 to 1 scale
    [train, pars] = stretch_cols_ind(train, -1, 1); 
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,unique(v)) = 0; % Remove NaNs from train set
    
    % Set train set to -1 to 1 scale
    [test2] = stretchWithGivenPars(test2, [-1 1], pars);
    test2(:,unique(v)) = 0; % Remove NaNs from test set
    
    fprintf('\nTraining SVM model for fold ' + string(f) + ' with ' + int2str(numK) + ' selected voxels...\n');
    pause(1);
    
    idx = find(ranks(:,f) <= numK); % Find indices that have a lower rank than numK
    train = train(:,idx);
    test2 = test2(:,idx);
    
    svm_model = svmtrain(gp, train, '-t 0 -c 1'); % -t 0 = linear SVM, -c 1 = cost value of 1
    svm_mod{f} = svm_model;
    
    % Test SVM model
    
    fprintf('\nPredicting conditions...\n');
    [svma_class(:,f), ~, ~] = svmpredict(CondClass', test2, svm_model);
    
    pause(1);
    svm_av(f) = length(find(svma_class(:,f) == CondClass')) ./ (nConditions);

end
% End cross-validation loop

% Proper output format
data{1} = train_set;
data{2} = test_set;
data{3} = anovas;

svmOut = [];
svmOut.models = svm_mod;    % SVM model for each cross-validation fold
svmOut.Aclass = svma_class; % Average classification
svmOut.pc = svm_av;         % Percentage correct single-block
svmOut.av = svm_av;         % Percentage correct average
svmOut.data = data;         % Parsed data used in classifier

% t test on performance across runs, within subject
[h,p,ci,stats] = ttest(svmOut.pc,1/nConditions,.05,'right');
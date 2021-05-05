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

% Parse for cross-validation cycles
[train_set, test_set, anovas] = parse_runs_surf(betasC);

% Output variables
svm_ws = cell(1,nRuns);
svm_class = zeros(nPerRun*(nConditions),nRuns);
svma_class = zeros(nConditions,nRuns);
svm_cm = zeros(nConditions,nConditions,nRuns);
svm_pc = zeros(nRuns,1);
svm_av = zeros(nRuns,1);
svm_mod = cell(1,nRuns);

% Define coding variable for training - assumes LOOCV (leaves one run out)
gp = []; k = 1; l = nPerRun*(nRuns - 1);
for ii = 1:nConditions
    gp(k:l,1) = ii;
    k = k + nPerRun*(nRuns - 1);
    l = l + nPerRun*(nRuns - 1);
end

% Define coding variable for testing
gp_test = []; k = 1; l = nPerRun;
for ii = 1:nConditions
    gp_test(k:l,1) = ii;
    k = k + nPerRun;
    l = l + nPerRun;
end

acc_avg = []; % Initialize accuracies
ranks = []; % Initialize ranks
nMax = size(train_set,2); % Maximum number of ranks

% Recursive Feature Elimination: Cross-validation cycles in every run
for r = 1:size(train_set,3)

    % Define train set and test set for this cycle of cross-validation
    train = train_set(:,:,r);  
    test2 = test_set(:,:,r);
    
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
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,unique(v)) = 0; % Remove NaNs from train set
    
    % Set train set to -1 to 1 scale
    [in] = stretchWithGivenPars(in, [-1 1], pars); % Test set for conditions
    in(:,v) = 0; % Remove NaNs from test set
    
    fprintf('\nTraining SVM model for run ' + string(r) + ' with RFE...\n');
    pause(2);
    
    [ftRank,~] = ftSel_SVMRFECBR(train,gp); % Compute RFE
    if r ~= size(train_set,3)
        acc = [];
        % Loop to compute accuracies with respect to the value of k
        for k = 1:nMax
            idx = find(ftRank <= k);
            train_k = train(:,idx);
            svm_model_k = svmtrain(gp, train_k, '-t 0 -c 1');
            [~, acc_k, ~] = svmpredict(CondClass', in, svm_model_k);
            acc = cat(1,acc,acc_k(1));
        end
        acc_avg = cat(2,acc_avg,acc);
    end
    ranks = cat(2,ranks,ftRank');
    
end

acc_avg = mean(acc_avg,2);
maxAcc = max(acc_avg);
numK = find(round(acc_avg,3) == round(maxAcc,3));
numK = numK(end);

% Main loop - Cross-validation cycles in every run
for r = 1:size(train_set,3)  

    % Define train set and test set for this cycle of cross-validation
    train = train_set(:,:,r);  
    test2 = test_set(:,:,r);
    
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
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,unique(v)) = 0; % Remove NaNs from train set
    
    % Set train set to -1 to 1 scale
    [test2] = stretchWithGivenPars(test2, [-1 1], pars);
    test2(:,unique(v)) = 0; % Remove NaNs from test set
    [in] = stretchWithGivenPars(in, [-1 1], pars); % Test set for conditions
    in(:,v) = 0; % Remove NaNs from test set
    
    fprintf('\nTraining SVM model for run ' + string(r) + ' with ' + int2str(numK) + ' selected voxels...\n');
    pause(2);
    
    idx = find(ranks(:,r) <= numK); % Find indices that have a lower rank than numK
    train = train(:,idx);
    test2 = test2(:,idx);
    
    svm_model = svmtrain(gp, train, '-t 0 -c 1'); % -t 0 = linear SVM, -c 1 = cost value of 1
    svm_mod{r} = svm_model;
    
    % Get the weights from model
    svm_weights = svm_DefineWeights(svm_model);  
    svm_ws{r} = svm_weights;
    
    % Test SVM model
    
    fprintf('\n\nPredicting trials/blocks...\n');
    pause(2);
    [svm_class(:,r), accuracy, dec] = svmpredict(gp_test, test2, svm_model);
    
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
    fprintf('\nPredicting conditions...\n');
    pause(2);
    [svma_class(:,r), accuracy, dec] = svmpredict(CondClass', in, svm_model);

    svm_av(r) = length(find(svma_class(:,r) == CondClass')) ./ (nConditions);

end
% End cross-validation loop

% Proper output format
data{1} = train_set;
data{2} = test_set;
data{3} = anovas;

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
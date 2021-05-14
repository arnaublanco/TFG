% Function that computes Support Vector Machine for the betas.
%  INPUT:
%   · betasC: Matrix containing the betas.
%   · CondClass: Stimuli (1: Forest, 2: People, 3: Traffic)
%   · permGP: 0 for normal analysis; 1 for randomization of labels
%  OUTPUT:
%   · knnOut: MATLAB object containing the results of the KNN.

function [knnOut] = singleKNN(betasC, CondClass, permGP)

nConditions = size(betasC,3);
nPerRun = size(betasC,1);
nRuns = size(betasC,4);

% Parse for cross-validation cycles
[train_set, test_set] = parse_runs_surf_blocks(betasC);

acc_avg = []; % Initialize average accuracy
% Cross-validation with KNN
for f = 1:size(train_set,3)  

    % Define train set and test set for this cycle of cross-validation
    train = train_set(:,:,f);  
    test2 = test_set(:,:,f);
    
    % Define coding variable for training - assumes LOOCV (leaves one run out)
    gp = []; k = 1; l = nPerRun*nRuns-1;
    for ii = 1:nConditions
        gp(k:l,1) = ii;
        k = k + nPerRun*nRuns-1;
        l = l + nPerRun*nRuns-1;
    end
    
    % Permutation of label vector - this will allow the creation of a randomization
    % distribution for good comparison purposes.
    if(permGP == 1)
        ff = randperm(length(gp)); % Create a random vector every time
        gp = gp(ff);
    end
    
    % If the training set and the labeling set do not have the dimensions,
    % then show an error.
    if(size(train,1) ~= size(gp,1))
        error('Training and gp vector mismatch');
    end
    
    % Normalize train and test sets
    [train, pars] = stretch_cols_ind(train, -1, 1); 
    [test2] = stretchWithGivenPars(test2, [-1 1], pars);
    
    % Remove NaNs in train and test sets
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,unique(v)) = []; % Remove NaNs from train set
    test2(:,unique(v)) = []; % Remove NaNs from test set
    
    acc = [];
    nMax = 100;
    for k = 1:nMax
        KNN_k = fitcknn(train,gp,'NumNeighbors',k,'Distance','Euclidean','Standardize',1);
        labels_k = predict(KNN_k,test2);
        acc_k = sum(CondClass' == labels_k)/size(labels_k,1);
        acc = cat(1,acc,acc_k);
    end
    acc_avg = cat(2,acc_avg,acc);
    
end
acc_avg = mean(acc_avg,2);
numK = find(round(acc_avg,3) == round(max(acc_avg),3));
numK = numK(end);

% Cross-validation with KNN
accuracy_c = [];
for f = 1:size(train_set,3)  

    % Define train set and test set for this cycle of cross-validation
    train = train_set(:,:,f);  
    test2 = test_set(:,:,f);
    
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
    
    % Permutation of label vector - this will allow the creation of a randomization
    % distribution for good comparison purposes.
    if(permGP == 1)
        ff = randperm(length(gp)); % Create a random vector every time
        gp = gp(ff);
    end
    
    % If the training set and the labeling set do not have the dimensions,
    % then show an error.
    if(size(train,1) ~= size(gp,1))
        error('Training and gp vector mismatch');
    end
    
    % Normalize train and test sets
    [train, pars] = stretch_cols_ind(train, -1, 1); 
    [test2] = stretchWithGivenPars(test2, [-1 1], pars);
    
    % Remove NaNs in train and test sets
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,unique(v)) = []; % Remove NaNs from train set
    test2(:,unique(v)) = []; % Remove NaNs from test set
    
    fprintf(['\nTraining SVM model for fold ',int2str(f)]);
    pause(1);
    
    % Fit KNN with train set
    KNN = fitcknn(train,gp,'NumNeighbors',numK,'Distance','Euclidean','Standardize',1);
    
    % Predict labels condition-wise
    fprintf('\nPredicting conditions...\n');
    labels_cond = predict(KNN,test2);
    
    % Compute accuracy condition-wise
    acc_cond = sum(CondClass' == labels_cond)/size(labels_cond,1);
    fprintf(['\nPrediction accuracy: ',num2str(acc_cond*100),'\n']);
    accuracy_c = cat(2,accuracy_c,acc_cond);
    
end

data{1} = train_set;
data{2} = test_set;

knnOut = [];
knnOut.pc = accuracy_c;
knnOut.av = accuracy_c;
knnOut.data = data;

end
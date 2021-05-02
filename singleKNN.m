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
[train_set, test_set, anovas] = parse_runs_surf(betasC);

acc_avg = []; % Initialize average accuracy
% Cross-validation with KNN
for r = 1:size(train_set,3)  

    % Define train set and test set for this cycle of cross-validation
    train = train_set(:,:,r);  
    test2 = test_set(:,:,r);
    
     % Define coding variable for training - assumes LOOCV (leaves one run out)
    gp = []; k = 1; l = nPerRun*(nRuns - 1);
    for ii = 1:nConditions
        gp(k:l,1) = ii;
        k = k + nPerRun*(nRuns - 1);
        l = l + nPerRun*(nRuns - 1);
    end
    
    % Permutation of label vector - this will allow the creation of a randomization
    % distribution for good comparison purposes.
    if(permGP == 1)
        f = randperm(length(gp)); % Create a random vector every time
        gp = gp(f);
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
    
    % Define coding variable for testing
    gp_test = []; k = 1; l = nPerRun;
    for ii = 1:nConditions
        gp_test(k:l,1) = ii;
        k = k + nPerRun;
        l = l + nPerRun;
    end
    
    % Normalize train and test sets
    [train, pars] = stretch_cols_ind(train, -1, 1); 
    [test2] = stretchWithGivenPars(test2, [-1 1], pars);
    [in] = stretchWithGivenPars(in, [-1 1], pars);
    
    % Remove NaNs in train and test sets
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,unique(v)) = []; % Remove NaNs from train set
    test2(:,unique(v)) = []; % Remove NaNs from test set
    in(:,unique(v)) = []; % Remove NaNs from test set
    
    acc = [];
    nMax = 100;
    for k = 1:nMax
        KNN_k = fitcknn(train,gp,'NumNeighbors',k,'Distance','Euclidean','Standardize',1);
        %labels_k = predict(KNN_k,test2);
        try predict(KNN_k,in)
            labels_k = predict(KNN_k,in);
        catch
            error('not working');
        end
        %acc_k = sum(gp_test == labels_k)/size(labels_k,1);
        acc_k = sum(CondClass' == labels_k)/size(labels_k,1);
        acc = cat(1,acc,acc_k);
    end
    acc_avg = cat(2,acc_avg,acc);
    
end
acc_avg = mean(acc_avg,2);
numK = find(round(acc_avg,3) == round(max(acc_avg),3));
numK = numK(1);

% Cross-validation with KNN
accuracy_b = []; accuracy_c = [];
for r = 1:size(train_set,3)  

    % Define train set and test set for this cycle of cross-validation
    train = train_set(:,:,r);  
    test2 = test_set(:,:,r);
    
     % Define coding variable for training - assumes LOOCV (leaves one run out)
    gp = []; k = 1; l = nPerRun*(nRuns - 1);
    for ii = 1:nConditions
        gp(k:l,1) = ii;
        k = k + nPerRun*(nRuns - 1);
        l = l + nPerRun*(nRuns - 1);
    end
    
    % Permutation of label vector - this will allow the creation of a randomization
    % distribution for good comparison purposes.
    if(permGP == 1)
        f = randperm(length(gp)); % Create a random vector every time
        gp = gp(f);
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
    
    % Define coding variable for testing
    gp_test = []; k = 1; l = nPerRun;
    for ii = 1:nConditions
        gp_test(k:l,1) = ii;
        k = k + nPerRun;
        l = l + nPerRun;
    end
    
    % Normalize train and test sets
    [train, pars] = stretch_cols_ind(train, -1, 1); 
    [test2] = stretchWithGivenPars(test2, [-1 1], pars);
    [in] = stretchWithGivenPars(in, [-1 1], pars);
    
    % Remove NaNs in train and test sets
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,unique(v)) = []; % Remove NaNs from train set
    test2(:,unique(v)) = []; % Remove NaNs from test set
    in(:,unique(v)) = []; % Remove NaNs from test set
    
    % Fit KNN with train set
    KNN = fitcknn(train,gp,'NumNeighbors',numK,'Distance','Euclidean','Standardize',1);
    
    % Predict labels block-wise and condition-wise
    labels_block = predict(KNN,test2);
    labels_cond = predict(KNN,in);
    
    % Compute accuracy block-wise and condition-wise
    acc_block = sum(gp_test == labels_block)/size(labels_block,1);
    acc_cond = sum(CondClass' == labels_cond)/size(labels_cond,1);
    
    accuracy_b = cat(2,accuracy_b,acc_block);
    accuracy_c = cat(2,accuracy_c,acc_cond);
    
end

knn_pc = accuracy_b;
knn_av = accuracy_c;

data{1} = train_set;
data{2} = test_set;
data{3} = anovas;

knnOut = [];
knnOut.pc = knn_pc;
knnOut.av = knn_av;
knnOut.data = data;

end
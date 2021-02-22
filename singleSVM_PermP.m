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

% to do single SVM analysis - once for whole of POI - no-subsampling
% only use SVM since LDA will not compute (cov will be singular)
% thus these results show --- method not dependent on voxel sub-sampling
% for significance!!!!! 

addpath('/Users/blancoarnau/Documents/MATLAB/libsvm/matlab/') % Import SVM toolbox 

nTrials = p(1);
nConditions = p(2);
nPerRun = p(3);
nVert = size(train_set,2);  % Always take updated value from betasC matrix
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
    
    %% *** Single-block SVM prediction ***
    
    % Set training set to -1 to 1 scale
    [train, pars] = stretch_cols_ind(train, -1, 1); 
    
    % Train SVM model
    [~,v] = find(isnan(train)); % Find NaNs in train set
    train(:,v) = []; % Remove NaNs from train set
    svm_model = svmtrain(gp, train, '-t 0 -c 1'); % -t 0 = linear SVM, -c 1 = cost value of 1
    svm_mod{r} = svm_model;
    
    % Get the weights from model
    svm_weights = svm_DefineWeights(svm_model);  
    svm_ws{r} = svm_weights;

    % Define test coding variable
    gp_test = []; k = 1; l = nPerRun;
    for ii = 1:nConditions
        gp_test(k:l,1) = ii;
        k = k + nPerRun;
        l = l + nPerRun;
    end

    %% ** Test SVM model **
    % Set test set to -1 to 1 scale
    [test2] = stretchWithGivenPars(test2, [-1 1], pars);
    % [test2, pars] = stretch_cols_ind(test2, -1, 1);
    
    [~,v] = find(isnan(test2)); % Find NaNs in test set
    test2(:,v) = []; % Remove NaNs from test set
    [svm_class(:,r), accuracy, dec] = svmpredict(gp_test, test2, svm_model);

    % Compute performance on test runs
    f = zeros(nConditions);  % Confusion matrix 

    k = 1; l = nPerRun;
    for i = 1:nConditions
        for j = 1:nConditions
            f(i,j) = length(find(svm_class(k:l,r)==j));
        end
        k = k + nPerRun; l = l + nPerRun;
    end

    svm_cm(:,:,r) = f;
    svm_pc(r) = trace(f)/(nPerRun*(nConditions));  

    %% ** Test SVM Prediction of Condition Average in Test Run **
    
    [in] = stretchWithGivenPars(in, [-1 1], pars);  % Set on same scale
    CondClass = 1:3;
    [svma_class(:,r), accuracy, dec] = svmpredict(CondClass', in, svm_model);
    accuracy(isnan(accuracy)) = 0; % If accuracy is NaN, set to 0
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
svmOut.pc = svm_pc;         % Percentage correct single trial/block
svmOut.av = svm_av;         % Percentage correct average
svmOut.data = data;         % Parsed data used in classifier

% t-test on performance across runs, within subject
[h,p,ci,stats]= ttest(svmOut.pc,1/nConditions,.05,'right');

% t=[];
% p=[];
% 
% if(marker==1 || marker==2)
%         
%     % single block
%     m = mean(svm_pc);  %% one mean per subject
%     t = (mean(m)-1/3)./(std(m)./sqrt(length(m)));
% 
%     df = length(m)-1;
%     p = 1-tcdf(t,df);  %% one tailed p for m>1/3
% 
%     % average
%     m = mean(svm_av);  %% one mean per subject
%     t(2) = (mean(m)-1/3)./(std(m)./sqrt(length(m)));
% 
%     df = length(m)-1;
%     p(2) = 1 - tcdf(t(2),df);
%     
% elseif(marker == 3)
%     
%     % single block
%     for r=1:2
%         m=svm_pc(r,:);  %% one mean per subject
%         t(r)=(mean(m)-1/3)./(std(m)./sqrt(length(m)));
% 
%         df=length(m)-1;
%         p(r)=1-tcdf(t(r),df);
%     end
% 
%     % average
%     for r=1:2
%         m=(svm_av(r,:));  %% one mean per subject
%         tAv(r)=(mean(m)-1/3)./(std(m)./sqrt(length(m)));
% 
%         df=length(m)-1;
%         pAv(r)=1-tcdf(tAv(r),df);
%     end
% end
%     
% 
% 
% outname=sprintf('SingleSVM_BlockDesign_Patch%d_Marker%d_zBetas%d_notUni%d.mat', Patch_ind,marker,zBetas,notUni);
% 
% save(outname);


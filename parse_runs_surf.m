% Function that returns the train set, test test and the ANOVA tables of
% a given matrix of betas.
%  INPUT:
%   · mask_vox2: Matrix of betas
%  OUTPUT:
%   · train_set: Training set (last iteration in CV)
%   · test_set: Test set (last iteration in CV)
%   · anovas: ANOVA tables (one per run) of the training set.

function [train_set, test_set, anovas] = parse_runs_surf(mask_vox2)

[nPerRun, nVert, nConditions, nRuns] = size(mask_vox2);

% The function 'nchoosek' will return a matrix with all possible
% combinations with 'nRuns' choose 'nRuns' - 1. Every number represents what run
% is taken out for training.
code = nchoosek(1:nRuns,nRuns-1); % Train on n-1 elements, and test the nth (LOOCV)

count = zeros(nRuns); % Count the number of times a given run appears in the permutation
train_set = zeros((nRuns-1)*nPerRun*nConditions, nVert, nRuns);
test_set = zeros(nConditions*nPerRun, nVert, nRuns);

% Leave-one-out cross-validation
for r = 1:nRuns
    
    training = []; % Collapse across runs
    for i = 1:nRuns-1
        training = cat(1,training,mask_vox2(:,:,:,code(r,i))); % Concatenate betas at permutated run
        count(r,code(r,i)) = count(r,code(r,i))+1; % Count that permutated run
    end

    left = find(count(r,:) == 0);
    t = mask_vox2(:,:,:,left);  % Find non-used run (i.e. the one left)
    
    % Collapse over conditions for both training and test %
    
    % Training
    train = [];
    for i = 1:nConditions
        train = cat(1,train,training(:,:,i));
    end

    % Testing
    test = [];
    for i = 1:nConditions
        test = cat(1,test,t(:,:,i));
    end
    
    % Store cross-validation fold
    train_set(:,:,r) = train;
    test_set(:,:,r) = test;
    
    % Add ANOVA for each training set
    gp = []; k = 1; l = nPerRun*(nRuns-1);
    for ii = 1:nConditions
        gp(k:l,1) = ii;
        k = k+nPerRun*(nRuns-1);
        l = l+nPerRun*(nRuns-1);
    end
    
    [anovas(:,:,r)] = voi_ANOVA(squeeze(train_set(:,:,r)), gp); % Compute ANOVAs per run of training set 
    
end



    
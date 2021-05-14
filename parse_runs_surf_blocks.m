% Function that returns the train set, test test and the ANOVA tables of
% a given matrix of betas.
%  INPUT:
%   · mask_vox2: Matrix of betas
%  OUTPUT:
%   · train_set: Training set (last iteration in CV)
%   · test_set: Test set (last iteration in CV)
%   · anovas: ANOVA tables (one per run) of the training set.

function [train_set, test_set] = parse_runs_surf_blocks(mask_vox2)

[nPerRun, nVert, nConditions, nRuns] = size(mask_vox2);

% The function 'nchoosek' will return a matrix with all possible
% combinations with nRuns choose 'nRuns' - 1. Every number represents what run
% is taken out for training.

nFolds = nPerRun*nRuns;

code = nchoosek(1:nFolds,nFolds-1); % Train on n-1 elements, and test the nth (LOOCV)

% Initialization of train set and test set
train_set = zeros((nFolds - 1)*nConditions, nVert, nFolds);
test_set = zeros(nConditions, nVert, nFolds);

% Leave-one-out cross-validation
for f = 1:nFolds
    
    % Collapse across runs
    collapse_runs = [];
    for i = 1:nRuns
        collapse_runs = cat(1,collapse_runs,mask_vox2(:,:,:,i)); % Concatenate betas at permutated fold
    end

    % Collapse over conditions
    collapse_cond = [];
    for i = 1:nConditions
        collapse_cond = cat(1,collapse_cond,collapse_runs(:,:,i));
    end
    
    % Select blocks for the test set (3 blocks in total)
    test = []; count = []; left = setdiff(1:nFolds,code(f,:));
    for c = 1:nConditions
       count = cat(1,count,(c-1)*nPerRun*nRuns+left);
       test = cat(1,test,collapse_cond((c-1)*nPerRun*nRuns+left,:));
    end
    
    % Select the rest of the blocks for the train set
    train = collapse_cond(setdiff(1:nPerRun*nRuns*nConditions,count),:);
    
    % Compute labels gp for ANOVA
    gp = []; k = 1; l = nPerRun*nRuns-1;
    for ii = 1:nConditions
        gp(k:l,1) = ii;
        k = k+nPerRun*nRuns-1;
        l = l+nPerRun*nRuns-1;
    end
    
    % Store cross-validation fold
    train_set(:,:,f) = train; 
    test_set(:,:,f) = test;
    
    %[anovas(:,:,f)] = voi_ANOVA(squeeze(train_set(:,:,f)), gp); % Compute ANOVAs per block of training set 
    
end



    
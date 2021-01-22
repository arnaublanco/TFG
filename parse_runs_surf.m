% Function that returns the train set, test test and the ANOVA tables of
% a given matrix of betas.
%  INPUT:
%   · mask_vox2: Matrix of betas
%  OUTPUT:
%   · train_set: Training set
%   · test_set: Test set
%   · anovas: ANOVA tables (one per run) of the training set.

function [train_set, test_set, anovas] = parse_runs_surf(mask_vox2)

[nPerRun, nVert, nConditions, nRuns] = size(mask_vox2);

% The function ´nchoosek´ will return a matrix where each row
% contains numbers from 1 to ´nRuns´. Every number represents what run
% is taken out for training.
code = nchoosek(1:nRuns,nRuns-1);  % Train on n-1 elements, and test the nth (LOOCV)

count = zeros(nRuns);
train_set = zeros((nRuns-1)*nPerRun*nConditions, nVert, nRuns);
test_set = zeros(nConditions*nPerRun, nVert, nRuns);

% Leave-one-out cross-validation
for r = 1:nRuns
    
    tmp = zeros(nPerRun,nVert,nConditions); % Collapse across runs
    
    for i = 1:nRuns-1
        tmp = cat(1,tmp,mask_vox2(:,:,:,code(r,i)));
        count(r,code(r,i)) = count(r,code(r,i))+1;
    end

    training = tmp(nPerRun+1:end,:,:);  %% should be nPerRun*(nRuns-1) by nVert by nCOnditions

    t = find(count(r,:) == 0);
    test1 = mask_vox2(:,:,:,t);  % Find non-used run
    
    % Collapse over conditions for both training and test %
    
    % Training
    tmp1 = zeros(nPerRun*(nRuns-1),nVert);
    for i = 1:nConditions
        tmp1 = cat(1,tmp1,training(:,:,i));
    end
    
    train = tmp1(nPerRun*(nRuns-1)+1:end,:);  

    % Testing
    tmp2 = zeros(nPerRun,nVert);
    for i = 1:nConditions
        tmp2 = cat(1,tmp2,test1(:,:,i));
    end
    
    test2 = tmp2(nPerRun+1:end,:);
    
    train_set(:,:,r) = train;
    test_set(:,:,r) = test2;
    
    % Add ANOVA for each training set (i.e. not on all data, just training)
    gp = []; k = 1; l = nPerRun*(nRuns-1);
    for ii = 1:nConditions
        gp(k:l,1) = ii;
        k = k+nPerRun*(nRuns-1);
        l = l+nPerRun*(nRuns-1);
    end
    
    [anovas(:,:,r)] = voi_ANOVA(squeeze(train_set(:,:,r)), gp); % Compute ANOVAs per run of training set 
    
end



    
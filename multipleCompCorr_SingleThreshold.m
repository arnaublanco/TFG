% Compute p-values

function multipleCompCorr_SingleThreshold(classifier)

if classifier == 1
    folder = 'SVM';
elseif classifier == 2
    folder = 'SVM_RFE';
else
    folder = 'KNN';
end
load(['Results/',folder,'/groupResults_Perm.mat']); % Load group results with permutation
nPerm = 1000;

% Accuracy for each ROI per single block/trial
allApc1 = meanApc.Auditory;
allApc1(:,2) = meanApc.Motor;
allApc1(:,3) = meanApc.V1;
allApc1(:,4) = meanApc.V2;

% Accuracy for each ROI per condition
allSpc1 = meanSpc.Auditory;
allSpc1(:,2) = meanSpc.Motor;
allSpc1(:,3) = meanSpc.V1;
allSpc1(:,4) = meanSpc.V2;

% Find maximum accuracy in both arrays - maximum accuracy distribution
maxApc1 = max(allApc1,[],2); 
maxSpc1 = max(allSpc1,[],2); 

% Calculate p-value for each ROI using the maximum accuracy distribution
pPermCorr_Apc.Auditory = length(find(maxApc1 >= meanObsApc.Auditory)) ./nPerm;
pPermCorr_Apc.Motor = length(find(maxApc1 >= meanObsApc.Motor)) ./nPerm;
pPermCorr_Apc.V1 = length(find(maxApc1 >= meanObsApc.V1)) ./nPerm;
pPermCorr_Apc.V2 = length(find(maxApc1 >= meanObsApc.V2)) ./nPerm;

pPermCorr_Spc.Auditory = length(find(maxSpc1 >= meanObsSpc.Auditory)) ./nPerm;
pPermCorr_Spc.Motor = length(find(maxSpc1 >= meanObsSpc.Motor)) ./nPerm;
pPermCorr_Spc.V1 = length(find(maxSpc1 >= meanObsSpc.V1)) ./nPerm;
pPermCorr_Spc.V2 = length(find(maxSpc1 >= meanObsSpc.V2)) ./nPerm;


save(['Results/',folder,'/groupResults_SingleThreshold_corrected'],'pPermCorr_Apc','pPermCorr_Spc');


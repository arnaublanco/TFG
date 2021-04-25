% This is the later version incorporating the non-biased computation of the confusion matrices
function run_single_sub_CollapseHem_get_unbiased_conf_mats(subject,Patch_ind,CondClass,POIfile_ind)

permGP = 0; % 0: Standard analysis; 1: Permutation analysis

% Prepare data - get timeseries, do GLM etc
[outDL] = readDataMtcPoi(subject, Patch_ind, 'LH', CondClass, POIfile_ind);
[outDR] = readDataMtcPoi(subject, Patch_ind, 'RH', CondClass, POIfile_ind);

% Concatenate across vertice dimension - hemispheres
outD = [];
f1 = size(outDL.betasC);
f2 = size(outDR.betasC);
if(f1(1) == f2(1) && f1(3) == f2(3) && f1(4) == f2(4))
    betasC = cat(2,outDL.betasC,outDR.betasC);
    betas = cat(2,outDL.betas,outDR.betas);
    S{3} = outDL.S{3};
    S{3}(4) = f1(2)+f2(2);
    f = []; 
    f.CondClass = CondClass;
    S{5} = f;
    outD.betasC = betasC;
    outD.betas = betas;
    outD.S = S;
    
else
    error('Hemispheric Data mismatch');
end

% Perform classification
[svmOut] = singleSVM(outD, permGP);

% t test across runs
nConditions=length(outD.S{5}.CondClass);  %% conditions being classified
[hST,pST,ci,statsST]=ttest(svmOut.pc,1/nConditions,.05,'right'); % single trial
[hAV,pAV,ci,statsAV]=ttest(svmOut.av,1/nConditions,.05,'right'); % average

% to save the output
outname = sprintf('%s_UnbiasedConfMatsCollapseHem_Patch%d_POI%d.mat', subject, Patch_ind, POIfile_ind);
save(outname, 'subject','Patch_ind','permGP','outD','svmOut');

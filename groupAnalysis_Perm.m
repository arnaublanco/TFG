% Obtain the results of the permutation analysis across subjects

function groupAnalysis_Perm(classifier)

subjects = 8;
patches = 2;
nPerm = 100;

if classifier == 1
    folder = 'SVM';
elseif classifier == 2
    folder = 'SVM_RFE';
else
    folder = 'KNN';
end

for subject_nr = 1:length(subjects)
    if subject_nr ~= 6
        for poifile = 1:3
            % Visual areas
            if poifile == 3
                    pooledApc_Poi3 = repmat(struct(),[1 patches]);
                    pooledSpc_Poi3 = repmat(struct(),[1 patches]);
                for patch = 1:patches
                    filename = ['Results/',folder,'/sub-0',int2str(subject_nr),'_MainAnalysis_CollapseHem_Patch',num2str(patch),'_POI',num2str(poifile)];
                    load(filename);
                    pooledApc_Poi3(patch).permAcc(1:nPerm,subject_nr) = Apc;
                    pooledSpc_Poi3(patch).permAcc(1:nPerm,subject_nr) = Spc;

                    pooledApc_Poi3(patch).obsAcc(1,subject_nr) = mean(OutObs.av);
                    pooledSpc_Poi3(patch).obsAcc(1,subject_nr) = mean(OutObs.pc);   
                end
            % Motor cortex
            elseif poifile == 2
                patch = 1;
                pooledApc_Poi2 = struct();
                pooledSpc_Poi2 = struct();
                filename = ['Results/',folder,'/sub-0',int2str(subject_nr),'_MainAnalysis_CollapseHem_Patch',num2str(patch),'_POI',num2str(poifile)];
                load(filename);
                
                pooledApc_Poi2.permAcc(1:nPerm,subject_nr) = Apc;
                pooledSpc_Poi2.permAcc(1:nPerm,subject_nr) = Spc;

                pooledApc_Poi2.obsAcc(1,subject_nr) = mean(OutObs.av);
                pooledSpc_Poi2.obsAcc(1,subject_nr) = mean(OutObs.pc);
            % Auditory cortex
            else
                patch = 1;
                pooledApc_Poi1 = struct();
                pooledSpc_Poi1 = struct();
                filename = ['Results/',folder,'/sub-0',int2str(subject_nr),'_MainAnalysis_CollapseHem_Patch',num2str(patch),'_POI',num2str(poifile)];
                load(filename);
                pooledApc_Poi1.permAcc(1:nPerm,subject_nr) = Apc;
                pooledSpc_Poi1.permAcc(1:nPerm,subject_nr) = Spc;
                % For the observed values
                pooledApc_Poi1.obsAcc(1,subject_nr) = mean(OutObs.av);
                pooledSpc_Poi1.obsAcc(1,subject_nr) = mean(OutObs.pc);

            end
        end
    end
end

save(['Results/',folder,'/pooledPermResults.mat'],'pooledApc_Poi3','pooledSpc_Poi3','pooledApc_Poi2','pooledSpc_Poi2','pooledApc_Poi1','pooledSpc_Poi1');

meanApc = struct();
meanSpc = struct();
meanObsApc = struct();

meanApc.Auditory = mean(pooledApc_Poi1(1).permAcc,2);
meanApc.Motor = mean(pooledApc_Poi2(1).permAcc,2);
meanApc.V1 = mean(pooledApc_Poi3(1).permAcc,2);
meanApc.V2 = mean(pooledApc_Poi3(2).permAcc,2);

meanSpc.Auditory = mean(pooledSpc_Poi1(1).permAcc,2);
meanSpc.Motor = mean(pooledSpc_Poi2(1).permAcc,2);
meanSpc.V1 = mean(pooledSpc_Poi3(1).permAcc,2);
meanSpc.V2 = mean(pooledSpc_Poi3(2).permAcc,2);

meanObsApc.Auditory = mean(pooledApc_Poi1(1).obsAcc,2);
meanObsApc.Motor = mean(pooledApc_Poi2(1).obsAcc,2);
meanObsApc.V1 = mean(pooledApc_Poi3(1).obsAcc,2);
meanObsApc.V2 = mean(pooledApc_Poi3(2).obsAcc,2);

meanObsSpc.Auditory = mean(pooledSpc_Poi1(1).obsAcc,2);
meanObsSpc.Motor = mean(pooledSpc_Poi2(1).obsAcc,2);
meanObsSpc.V1 = mean(pooledSpc_Poi3(1).obsAcc,2);
meanObsSpc.V2 = mean(pooledSpc_Poi3(2).obsAcc,2);

group_pPerm_Apc.Auditory = length(find(meanApc.Auditory >= meanObsApc.Auditory)) ./ nPerm;
group_pPerm_Apc.Motor = length(find(meanApc.Motor >= meanObsApc.Motor)) ./ nPerm;
group_pPerm_Apc.V1 = length(find(meanApc.V1 >= meanObsApc.V1)) ./ nPerm;
group_pPerm_Apc.V2 = length(find(meanApc.V2 >= meanObsApc.V2)) ./ nPerm;

group_pPerm_Spc.Auditory = length(find(meanSpc.Auditory >= meanObsSpc.Auditory)) ./ nPerm;
group_pPerm_Spc.Motor = length(find(meanSpc.Motor >= meanObsSpc.Motor)) ./ nPerm;
group_pPerm_Spc.V1 = length(find(meanSpc.V1 >= meanObsSpc.V1)) ./ nPerm;
group_pPerm_Spc.V2 = length(find(meanSpc.V2 >= meanObsSpc.V2)) ./ nPerm;

save(['Results/',folder,'/groupResults_Perm.mat'],'meanSpc','meanApc','meanObsApc','meanObsSpc','group_pPerm_Apc','group_pPerm_Spc');





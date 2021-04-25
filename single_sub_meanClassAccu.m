% Computes the mean classification accuracy results for each subject
% across cross-validation cycles.
function single_sub_meanClassAccu

for subject = 1:length(subjects)
    for poi = 1:3 
        % Visual areas
        if poi == 3
            for patch = 1:2

                filename = ['.Results/sub-0',int2str(subject),'_MainAnalysis_CollapseHem_Patch',int2str(patch),'_POI',int2str(poi),'.mat'];
                load(filename);
                meanAccu_singleblock = mean(svmOutObs.pc);
                meanAccu_average = mean(svmOutObs.av);
                semAccu_singleblock = std(svmOutObs.pc)/sqrt(length(svmOutObs.pc));
                semAccu_average = std(svmOutObs.av)/sqrt(length(svmOutObs.pc));           
                results_visual(patch,:) = [meanAccu_singleblock,meanAccu_average, semAccu_singleblock, semAccu_average];
            end
            save(['.Results/sub-0',int2str(subject),'meanClassAccu_visual_cortex'],'results_visual');
        % Motor cortex
        elseif poi == 2
            patch = 1;

            filename = ['.Results/sub-0',int2str(subject),'_MainAnalysis_CollapseHem_Patch',int2str(patch),'_POI',int2str(poi),'.mat'];
            load(filename);           
            meanAccu_singleblock = mean(svmOutObs.pc);
            meanAccu_average = mean(svmOutObs.av);
            semAccu_singleblock = std(svmOutObs.pc)/sqrt(length(svmOutObs.pc));
            semAccu_average = std(svmOutObs.av)/sqrt(length(svmOutObs.pc));

            results_motor(patch,:) = [meanAccu_singleblock,meanAccu_average, semAccu_singleblock, semAccu_average];
            save(['.Results/sub-0',int2str(subject),'meanClassAccu_motor_cortex'],'results_motor');
        % Auditory cortex
        else
           patch = 1;

           filename = ['.Results/sub-0',int2str(subject),'_MainAnalysis_CollapseHem_Patch',int2str(patch),'_POI',int2str(poi),'.mat'];
           load(filename);
           meanAccu_singleblock = mean(svmOutObs.pc);
           meanAccu_average = mean(svmOutObs.av);
           semAccu_singleblock = std(svmOutObs.pc)/sqrt(length(svmOutObs.pc));
           semAccu_average = std(svmOutObs.av)/sqrt(length(svmOutObs.pc));

           results_auditory(patch,:) = [meanAccu_singleblock,meanAccu_average, semAccu_singleblock, semAccu_average];
           save(['.Results/sub-0',int2str(subject),'meanClassAccu_auditory_cortex'],'results_auditory');
        end
    end
end






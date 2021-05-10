% Function that computes the mean classification accuracy results for each subject
% across cross-validation cycles.
%  INPUT:
%   · classifier: Type of classifier.

function single_sub_meanClassAccu(classifier)

subjects = 8;
patches = 2;

if classifier == 1
    folder = 'SVM';
elseif classifier == 2
    folder = 'SVM_RFE';
else
    folder = 'KNN';
end

for subject = 1:length(subjects)
    if subject ~= 6
        for poi = 1:3 
            % Visual areas
            if poi == 3
                results_visual = zeros(patches,4);
                for patch = 1:patches

                    filename = ['Results/',folder,'/sub-0',int2str(subject),'_MainAnalysis_CollapseHem_Patch',int2str(patch),'_POI',int2str(poi),'.mat'];
                    load(filename);
                    meanAccu_singleblock = mean(OutObs.pc);
                    meanAccu_average = mean(OutObs.av);
                    semAccu_singleblock = std(OutObs.pc)/sqrt(length(OutObs.pc));
                    semAccu_average = std(OutObs.av)/sqrt(length(OutObs.pc));           
                    results_visual(patch,:) = [meanAccu_singleblock, meanAccu_average, semAccu_singleblock, semAccu_average];
                end
                save(['Results/',folder,'/sub-0',int2str(subject),'meanClassAccu_visual_cortex'],'results_visual');
            % Motor cortex
            elseif poi == 2
                patch = 1;

                filename = ['Results/',folder,'/sub-0',int2str(subject),'_MainAnalysis_CollapseHem_Patch',int2str(patch),'_POI',int2str(poi),'.mat'];
                load(filename);           
                meanAccu_singleblock = mean(OutObs.pc);
                meanAccu_average = mean(OutObs.av);
                semAccu_singleblock = std(OutObs.pc)/sqrt(length(OutObs.pc));
                semAccu_average = std(OutObs.av)/sqrt(length(OutObs.pc));

                results_motor = [meanAccu_singleblock, meanAccu_average, semAccu_singleblock, semAccu_average];
                save(['Results/',folder,'/sub-0',int2str(subject),'meanClassAccu_motor_cortex'],'results_motor');
            % Auditory cortex
            else
               patch = 1;

               filename = ['Results/',folder,'/sub-0',int2str(subject),'_MainAnalysis_CollapseHem_Patch',int2str(patch),'_POI',int2str(poi),'.mat'];
               load(filename);
               meanAccu_singleblock = mean(OutObs.pc);
               meanAccu_average = mean(OutObs.av);
               semAccu_singleblock = std(OutObs.pc)/sqrt(length(OutObs.pc));
               semAccu_average = std(OutObs.av)/sqrt(length(OutObs.pc));

               results_auditory = [meanAccu_singleblock, meanAccu_average, semAccu_singleblock, semAccu_average];
               save(['Results/',folder,'/sub-0',int2str(subject),'_meanClassAccu_auditory_cortex'],'results_auditory');
            end
        end
    end
end






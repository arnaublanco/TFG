% Function that computes classification for all the subjects separately.
% INPUT:
%   · type: 1 for volume-based data and 2 for surface-based data.
%   · classifier: 1 for SVM and 2 for SVM-RFE
%   · cv: 1 for cross-validation per runs and 2 for CV per blocks

function run_all_sub_main_analysis(type, classifier, cv)

    subjects = 8; % Number of subjects
    CondClass = [{'forest'}, {'people'}, {'traffic'}]; % Conditions to classify (stimuli): Forest, People and Traffic (1, 2 and 3)
    patches = 2; % Number of patches for the EVC

    for subject = 1:subjects
        if subject ~= 6 % Skip the 6th subject (corrupted file)
            for POIfile_ind = 1:3  % 3 POIs/ROIs: Auditory, Motor and EVC
                if POIfile_ind == 3
                    for patch = 1:patches
                        run_single_sub_CollapseHem_main_analysis_parallel(subject, patch, CondClass, POIfile_ind, type, classifier, cv);
                    end
                else
                    patch = 1;
                    run_single_sub_CollapseHem_main_analysis_parallel(subject, patch, CondClass, POIfile_ind, type, classifier, cv);
                end
            end
        end
    end

end
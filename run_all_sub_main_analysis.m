% Function that computes classification for all the subjects separately.
% INPUT:
%   · type: 1 for volume-based data and 2 for surface-based data.
%   · classifier: 1 for SVM, 2 for SVM-RFE and 3 for KNN
function run_all_sub_main_analysis(type, classifier)

    subjects = 8; % Number of subjects
    CondClass = [{'forest'}, {'people'}, {'traffic'}]; % Conditions to classify (stimuli): Forest, People and Traffic (1, 2 and 3)
    if type == 1 % Volume-based
        patches = 3; % Number of patches of the EVC (V1, V2 and V3)
    else % Surface-based
        patches = 2; % Number of patches of the EVC (V1 and V2)
    end

    for subject = 1:subjects
        if subject ~= 6 % Skip the 6th subject (corrupted file)
            for POIfile_ind = 1:3  % 3 POIs/ROIs: Auditory, Motor and EVC
                if POIfile_ind == 3
                    for patch = 1:patches
                        run_single_sub_CollapseHem_main_analysis_parallel(subject, patch, CondClass, POIfile_ind, type, classifier);
                    end
                else
                    patch = 1;
                    run_single_sub_CollapseHem_main_analysis_parallel(subject, patch, CondClass, POIfile_ind, type, classifier);
                end
            end
        end
    end

end
% Run the single SVM for all subjects

subjects = 8;
CondClass = [{'forest'}, {'people'}, {'traffic'}]; % Conditions to classify: Forest, People and Traffic (1, 2 and 3)


for subject = 1:subjects
    for POIfile_ind = 1:3    % 3 POIs: Aud, Motor and EVC
        if POIfile_ind == 3 % EVC: EVC, V1, V2, V3, V1 fovea, V1 periphery, V1 far periphery, V2 fovea, etc. (total=16)
            for patch = 1:16
                run_single_sub_CollapseHem_main_analysis_parallel(subject,patch,CondClass,POIfile_ind);
            end
        else
            patch = 1;
            run_single_sub_CollapseHem_main_analysis_parallel(subject,patch,CondClass,POIfile_ind);
        end
    end
end

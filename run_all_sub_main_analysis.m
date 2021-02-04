% Run the single SVM for all subjests
close all;
clear all;

subjects = 8;
CondClass = [{'forest'}, {'people'}, {'traffic'}]; % Conditions to classify (stimuli): Forest, People and Traffic (1, 2 and 3)
visualize = 0; % File visualization as the code is running (1: 'on', 0: 'off')

for subject = 1:subjects
    for POIfile_ind = 1:3   % 3 POIs: Aud, Motor and EVC
        if POIfile_ind == 3 % EVC: EVC, V1, V2, V3, V1 fovea, V1 periphery, V1 far periphery, V2 fovea, etc. (total=16)
            for patch = 1:2
                run_single_sub_CollapseHem_main_analysis_parallel(subject,patch,CondClass,POIfile_ind,visualize);
            end
        else
            patch = 1;
            run_single_sub_CollapseHem_main_analysis_parallel(subject,patch,CondClass,POIfile_ind,visualize);
        end
    end
end

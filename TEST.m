clear all;
close all;

%% Data acquistion

data_dir = '/Users/blancoarnau/Desktop/';
%_space-MNI152NLin2009cAsym_res-2_desc-preproc
%data_dir = 'C:/Users/Arnau/Documents/TFG/ds002715/';
nRuns = 1;
subject = 1;

figure;
for i=1:nRuns
    
    run = i;
    disp(run);

    if subject < 10
        s = "sub-0" + string(subject);
    else
        s = "sub-" + string(subject);
    end
    
    if subject < 2
        file = s + "_ses-mri_task-AVScenesshort_run-" + string(run) + "_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz";
    else
        file = s + "_ses-mri_task-AVSceneslong_run-" + string(run) + "_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz";
    end
    
    path = data_dir + s + "/" + file;

    V = niftiread(path);
    info = niftiinfo(path);

   %%% Data visualization

    x = 30;
    y = 30;
    z = 30;
    
    voxel = double(reshape(V(x,y,z,:),1,size(V,4)));
    plot(1:size(voxel,2),voxel);
    title('Voxel ('+string(x)+','+string(y)+','+string(z)+')');
    hold on

    %dm = ones(size(voxel,2),2);
    %dm(1:50,2) = 0;
    %[B,DEV,STATS] = glmfit(dm,voxel,'normal','constant','off');

    %y = B(1)*dm(:,1)' + B(2)*dm(:,2)';
    %plot(1:size(voxel,2),y);
    hold on

end

%% To open .tsv files

a = 'test.tsv';
readtsv(a);

%%

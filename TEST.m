nSubjects = 8;
nRuns = 4;
CondClass = [{'forest'}, {'people'}, {'traffic'}];
POIfile_ind = 1;

Vms = zeros(62,73,73,nSubjects*nRuns);
for s = 1:nSubjects
    for r = 1:4
        if s ~= 6
            disp('Subject '+string(s) + ' | Run ' + string(r));
            [input_file_info] = getFileInfo(s, '', CondClass, 1);
            main = niftiread(input_file_info.func_name{r});

    %         i = 1;
    %         V = mean(main,4);
    %         while(i <= size(V,3))
    %             Vm = V(:,:,i);
    %             imshow(Vm,[min(Vm,[],'all') max(Vm,[],'all')])
    %             title('y = '+string(i))
    %             hold on
    %             pause(0.05)
    %             if i == size(V,3)
    %                 i = 1;
    %             else
    %                 i = i + 1;
    %             end
    %         end

            V = mean(main,4);
            Vms(:,:,:,s+(r-1)) = V;
        end
    end
end

%% Visualize

figure;
i = 1;
pause(3);
while(i <= size(Vms,3))
    
    for n = 1:size(Vms,4)
        V1 = Vms(:,:,:,n);
        for m = 1:size(Vms,4)
            V2 = Vms(:,:,:,m);
            %Vdiff = abs(V2-V1);
            subplot(7,7,(n-1)*7+1+(m-1))
            %Vm = Vdiff(:,:,i);
            %imshow(Vm,[min(Vm,[],'all') max(Vm,[],'all')])
            title('V'+string(n)+' - V'+string(m));
            imagesc(V1(:,:,i));
            hold on
            im = imagesc(V2(:,:,i));
            if i == 1
                colorbar; caxis([0 2000]);
                axis off
            end
        end
    end
    
    pause(0.2)
    if i == size(Vms,3)
        i = 1;
    else
        i = i + 1;
    end
end

%%

figure;
axis off
pause(3);
for t = 1:100
    for i = 1:4
        Vi = Vms(:,:,30,i);
        for j = 1:4
            subplot(4,4,(i-1)*4+1+(j-1))
            if mod(t,2) == 0
                imshow(Vi,[min(Vi,[],'all') max(Vi,[],'all')]);
            else
                Vj = Vms(:,:,30,j);
                imshow(Vj,[min(Vj,[],'all') max(Vj,[],'all')]);
            end
            title('Run '+string(i)+'- Run '+string(j));
            hold on
        end
    end
    pause(0.5);
end
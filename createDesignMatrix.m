function createDesignMatrix(subj,nPreds,dir_file,dir_dm_b,dir_dm_c)

conds = [{'forest'} {'traffic'} {'people'}];
    
if subj == 1
    timepoints = 117;
else
    timepoints = 222;
end

C = readtsv(dir_file);

dm_cond = zeros(timepoints,size(conds,2)+1);
dm_block = zeros(timepoints,nPreds);
dm_cond(:,1) = 1;
dm_block(:,1) = 1;

counter = 2;
for i=1:size(conds,2)
   idx = find(strcmp([C{3}], conds{i}));

   for b = 1:size(idx,1)
      dm_cond(C{1}(idx(b)):C{1}(idx(b))+C{2}(idx(b)),i+1) = 1;
      dm_block(C{1}(idx(b)):C{1}(idx(b))+C{2}(idx(b)),counter) = 1; 
      counter = counter + 1;
   end

end

dlmwrite(dir_dm_c,dm_cond,'delimiter','\t');
dlmwrite(dir_dm_b,dm_block,'delimiter','\t');

end
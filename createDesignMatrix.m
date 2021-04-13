function createDesignMatrix(subj,nPreds,dir_file,dir_dm_b,dir_dm_c)

conds = [{'forest'} {'traffic'} {'people'}];
TR = 2;

C = readtsv(dir_file); 

% The first participant went through a shorter version of the experiment.
if subj == 1
    timepoints = 117;
    C{2}(end) = 12; % Correcting last cell, which should be 12
else
    timepoints = 222;
end

dm_cond = zeros(timepoints,size(conds,2)+1);
dm_block = zeros(timepoints,nPreds);
dm_cond(:,1) = 1;
dm_block(:,1) = 1;

t = 1:TR:timepoints*TR; % Time during the experiment
h = HRF(t); % Hemodynamic Response Function

counter = 2;
for i = 1:size(conds,2)
   idx = find(strcmp([C{3}], conds{i}));

   for b = 1:size(idx,1)
       
      begin = round(C{1}(idx(b))/TR);
      fin = round(C{2}(idx(b))/TR);
      
      dm_cond(begin+1:begin+fin,i+1) = 1;
      dm_block(begin+1:begin+fin,counter) = 1; 
      counter = counter + 1;
   end
   tmp = conv(dm_cond(:,i+1),h'); % Convolve ´dm_cond´ with HRF
   dm_cond(:,i+1) = tmp(1:timepoints); % Place the values in ´dm_cond´
end

% Convolve ´dm_block´ with HRF
for j = 2:size(dm_block,2)
   tmp = conv(dm_block(:,j),h');
   dm_block(:,j) = tmp(1:timepoints);
end

dlmwrite(dir_dm_c,dm_cond,'delimiter','\t');
dlmwrite(dir_dm_b,dm_block,'delimiter','\t');

end
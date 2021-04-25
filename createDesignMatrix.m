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

%counter = 2;
for i = 1:size(conds,2)
   idx = find(strcmp([C{3}], conds{i}));

   for b = 1:size(idx,1)
       
      begin = round(C{1}(idx(b))/TR);
      fin = round(C{2}(idx(b))/TR);
      
      dm_cond(begin+1:begin+fin,i+1) = 1;
   end
   tmp = conv(dm_cond(:,i+1),h'); % Convolve ´dm_cond´ with HRF
   dm_cond(:,i+1) = tmp(1:timepoints); % Place the values in ´dm_cond´
end

idx = strcmp(C{3},{'baseline'});
C{1}(idx) = [];
C{2}(idx) = [];
C{3}(idx) = [];

for i = 1:nPreds-1
    begin = round(C{1}(i)/TR);
    fin = round(C{2}(i)/TR);
    dm_block(begin+1:begin+fin,i+1) = 1;
    tmp = conv(dm_block(:,i+1),h');
    dm_block(:,i+1) = tmp(1:timepoints);
end

writematrix(dm_cond,dir_dm_c);
writematrix(dm_block,dir_dm_b);

end
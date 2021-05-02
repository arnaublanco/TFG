%%% This file is intended to be run separately from the main code.
% It generates the ROI files from the 'broadmann.nii.gz' file.
% Template in MNI space from https://www.nitrc.org/projects/mricron .

% Regions of Interest (ROIs):
% - Auditory cortex: Brodmann areas 41, 42 and 22.
% - Visual cortex: Brodmann areas 17, 18 and 19.
% - Motor cortex: Brodmann area 4.

dir = 'brodmann.nii.gz'; % NifTi file directory.

V = niftiread(dir); % Read NifTi file.

V_auditory = zeros(size(V));
V_motor = zeros(size(V));
V_visual = zeros(size(V));

% The ROI masks are binary matrices that contain 1s' where the region of
% interest is, and 0s' otherwise.
V_visual(V == 17) = 1;
V_visual(V == 18) = 2;
V_visual(V == 19) = 3;
V_auditory(V == 41 | V == 42 | V == 22) = 1;
V_motor(V == 4) = 1;

% Resize the masks so that the fit the NifTi files' resolution.
dim = [62 73 73]; % Dimensions of all the brain scans (width, height and depth)
V_visual = imresize3(V_visual,[dim(1) dim(2) dim(3)],'nearest');
V_auditory = imresize3(V_auditory,[dim(1) dim(2) dim(3)],'nearest');
V_motor = imresize3(V_motor,[dim(1) dim(2) dim(3)],'nearest');

% Save matrices as NifTi files.
niftiwrite(V_visual,'Visual.nii');
niftiwrite(V_auditory,'Auditory.nii');
niftiwrite(V_motor,'Motor.nii');
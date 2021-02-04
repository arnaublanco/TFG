close all;
clear all;

addpath(genpath('/Applications/freesurfer/matlab'));
path = '/Users/blancoarnau/Documents/freesurfer/sub-01/mri/brain.mgz';
c = MRIread(path);
volume = c.vol;
volshow(volume);
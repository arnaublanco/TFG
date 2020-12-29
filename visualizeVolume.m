% This function visualizes a volume, typically an fMRI image.
function visualizeVolume(V)
% INPUT:
%   - V: Volume to visualize.    

    % Normalization: values between 0 and 1 -> x 255 for gray-scale visualization
    vmax = max(V,[],'all'); % Maximum value in the whole volume
    vmin = min(V,[],'all'); % Minimum value in the whole volume

    Vm = ((V-vmin)./(vmax-vmin)).*255;
    volshow(Vm);
end


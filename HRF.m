function [hrf] = HRF(timePoints)
% Gamma pdf for the peak
peak_values = gampdf(timePoints, 6);
% Gamma pdf for the undershoot
undershoot_values = gampdf(timePoints, 12);
% Combine them
values = peak_values - 0.35 * undershoot_values;
% Scale max to 0.6
hrf = values / max(values) * 0.6;
end
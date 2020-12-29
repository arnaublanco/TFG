% Function that reads a .txt file and returns what it contains.
function [outFile] = readDM(path, nReg)
% INPUT:
%   - path: Path to DM file.
%   - nReg: Number of regressors (number of columns).

    fobj = fopen(path,'r'); % Open file
    file = textscan(fobj, repmat('%f',[1 nReg])); % Read formatted file
    outFile = zeros(size(file{1},1),nReg); % Initialize output file ´outFile´
    for i = 1:nReg
       outFile(:,i) = file{i}; % cell array -> double array
    end
    fclose(fobj); % Close file
end


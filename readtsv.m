% Function that reads .tsv files.
%  INPUT:
%   · path: Path to file.
%  OUTPUT:
%   · file: Retrieved file.

function [file] = readtsv(path)
    fobj = fopen(path);
    file = textscan(fobj, '%f %f %s', 'HeaderLines', 1);
    fclose(fobj);
end


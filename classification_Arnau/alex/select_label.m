function [lx, ly, iLab ] = select_label(x, y, label)
% Usage [lx,ly, iLab] = select_label(x, y, label)
%    x: dataset
%    y: labelset
%    label: label to extract
%       This _function returns the selected dataset corresponding to 
%       the labels defined, the labels themselves and the correspondin indexes
%
%    Alex Perera, eSense Systems 2002
if ((max(y(:,1))>1) && (size(y,2)==1)),
  y=yRow2yStd(y);
end

iLab=[];

%old method
%for il=label,
%  %iLab = [iLab; find(y(:,il) == 1 )];
%end

iLab = [ iLab; find(sum(y(:,label),2) ==1)];
ly = y(iLab,:);
lx = x(iLab,:);
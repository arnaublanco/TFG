function [lx, ly, iLab ] = selectdata(x, y, label)
% Usage [lx,ly, iLab] = selectdata(x, y, label)
%    x: dataset
%    y: labelset
%    label: label to extract
%           [ 1 2 3] or { 1, [2 3]}
%       This _function returns the selected dataset corresponding to 
%       the labels defined, the labels themselves and the correspondin indexes
%
%    Alex Perera, eSense Systems 2002

if ~iscell(label)
    iLab=[];
    for il=label,
        iLab = [iLab; find(y(:,il) == 1 )];
    end
    iLab=sort(iLab);
    ly = y(iLab,:);
    lx = x(iLab,:);
else
    u=0; tLab=[];
    nvecs = size(x,1);
    ty = zeros(nvecs,length(label));
    for ic=1:length(label),
        iLab = find(sum(y(:,label{ic}),2) >0);
        tLab = [tLab; iLab];
        ty(iLab,ic)=1;   
    end
    tLab=sort(tLab);
    ly = ty(tLab,:);
	if size(x,3) >1,
		lx = x(tLab,:,:);
	else
		lx = x(tLab,:);
	end
    iLab=tLab;
end
    
    
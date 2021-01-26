function datetq=knn1(data,yStd,k,xin,Moptions)
% K-nearest neighbours algorithm
%
% Usage: yclas = knn(x,y,k,xin)
%
% Inputs:
%   -data: dataset m x n matrix; m samples of n dimensions
%	-yStd: labels
%	-k: number of neighbours
%   -xin: sample to clasify
%
% Outputs:
%   -yclas: 

if nargin < 5
    Moptions=ispmetrics('options','Euclidean');
end

[npat, dim]=size(data);
[ntest, dimtest]=size(xin);
maxlab = size(yStd,2);

if (dim	~=  dimtest)
    disp(['La dimension del vector de entrada' ...
        'no coincide con la del vector de test' ]);
end


datetq=zeros(1,ntest);


dist=ispmetrics(data,xin,Moptions);

parfor nvecin=1:ntest,
    [vecord, indsort]=sort(dist(:,nvecin),'ascend');
    toty=sum([yStd(indsort(1:k),:); zeros(1,maxlab)]); %#ok<PFBNS>
    [yt,indexok]=max(toty);
    datetq(1,nvecin)=indexok;
end

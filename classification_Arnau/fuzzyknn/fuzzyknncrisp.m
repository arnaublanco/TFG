function [clas]=fuzzyknncrisp(data,yStd,xin,k)

% Usage: [clas]=fuzzyknncrisp(data,yStd,xin,k)
%
% Returns the maximum class of a fuzzy knn proccess.
% using m=2
%
% Also see: fuzzyknncrisp, fuzzyknnval
%
% March 2003 Rafa Rubio

m=2;


[npat dim]=size(data);

[ntest dimtest]=size(xin);

maxlab = size(yStd,2);



if (dim	~=  dimtest)
    disp('La dimension del vector de entrada no coincide con la del vector de test');
end

if (npat ~= size(yStd,1))
    disp(' La dimension del vector de entrad no coincide con la de yStd');    
end



for nvecin=1:ntest
    vec=xin(nvecin,:);
    for patac=1:npat,
        if sum(vec~=data(patac,:))
            d(patac)=sqrt(sum((vec-data(patac,:)).^2));
        else 
            d(patac)=realmin;
        end
    end
    indsort=[];
    [vecord, indsort]=sort(d);
    if k~=1
        dataStd(nvecin,:)=sum( diag(min(ones(1,k).*realmax,vecord(1:k).^-(2/(m-1)) )) * yStd(indsort(1:k),:) ) / (sum (min( ones(1,k).*realmax , vecord(1:k).^-(2/(m-1)) )));  
    else
        dataStd(nvecin,:)=diag(min(ones(1,k).*realmax,vecord(1:k).^-(2/(m-1)) )) * yStd(indsort(1:k),:)  / (sum (min( ones(1,k).*realmax , vecord(1:k).^-(2/(m-1)) )));  
    end
    %%%%%%%%%%%%%%%%%%%% Decision %%%%%%%%%%%%%%%%%%%%%%
    [vv,ii]=max(dataStd(nvecin,:));
    clas(nvecin)=ii;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end










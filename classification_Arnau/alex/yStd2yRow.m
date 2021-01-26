function ySt=yStd2yRow(yStd)
% Usage yRow=yStd2yRow(yStd)
%    This fucntion converts ISPG Standard (binary) labels style 
%    to Row label style
%     
%  Alex Perera   ISPG/eSenseSystems 
%  15 May 2001

ySt=[];
[iyStd ,b]=find(yStd>0.5);%voy a cambiar yStd<0.5
%iyStd
%ySt=ones(max(iyStd-1),1);
ySt(iyStd)=b;
ySt=ySt';
%[ycer,a]=find(ySt==0)


   
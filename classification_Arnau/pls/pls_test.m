import classification.pls.*;

clc;

Xt=load('prova_pls_Xt.txt');
Yt=load('prova_pls_Yt.txt');

model.Xt=Xt;
model.Yt=Yt;
model.npls=5;
model.prep=1;
model.cl=95;

mtrain=pls_train(model);
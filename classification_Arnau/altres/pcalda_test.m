
clc;

Xt=load('prova_pcalda_train.txt');
it=load('prova_pcalda_itrain.txt');
Xval=load('prova_pcalda_val.txt');

model.Xt=Xt;
model.npc=6;
model.it=it;
model.prep=1;
model.cl=95;

mtrain=pcalda_train(model);
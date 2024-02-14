close;clear;clc;

load ..\data\P14_DYN_50_P35.mat
load ..\data\P14model-ocv.mat

data.temp    = 35; % temperature
data.script1 = DYNData.script1;
data.script2 = DYNData.script2;
data.script3 = DYNData.script3;

model = processDynamic(data,model,1,1);
save('..\data\P14model_dynamic.mat','model')
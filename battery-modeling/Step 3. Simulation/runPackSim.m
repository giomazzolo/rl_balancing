% Run random pack simulaitons

Ns = 3; % Number of cells
Nc = 100; % Number of cycles to simulate

cycleFiles = {'nyccPower.mat','uddsPower.mat','us06Power.mat','hwfetPower.mat'}; % drive cycles power demand files
cycleInd = 2;

if ~isdeployed; addpath ..\data; end

load("..\data\P14model_dynamic.mat", "model"); % Data Obtained from Dynamic Processing

% Dynamic model with milti temperature feature and different cell models
% fname = 'A123modeldyn.json'; 
% fid = fopen(fname); 
% raw = fread(fid,inf); 
% str = char(raw'); 
% fclose(fid); 
% model = jsondecode(str);

% Enable/Disble randomized variables
tOpts = 0; % random cell teperatures
qOpt = 1; % random cell capacities
rOpt = 1; % random cell r0
sdOpt = 1; % random cell self discharge 
cOpt = 1; % random cell-coulombic-efficiency values
lOpt = 1; % random cell leakage currents
rOpts = [tOpts, qOpt, rOpt, sdOpt, cOpt, lOpt];


packData = simRandPack(Ns,Nc,cycleFiles{cycleInd},model,rOpts);


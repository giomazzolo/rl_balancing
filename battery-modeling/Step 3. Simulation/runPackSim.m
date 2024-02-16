% Run random pack simulaitons

Ns = 6; % Number of cells
Nc = 100; % Number of cycles to simulate

cycleFiles = {'nyccPower.mat','uddsPower.mat','us06Power.mat','hwfetPower.mat'}; % drive cycles power demand files
cycleInd = 2;

load("..\data\P14model_dynamic.mat", "model"); % Data Obtained from Dynamic Processing

% Enable/Disble randomized variables
tOpts = 0; % random cell teperatures
qOpt = 0; % random cell capacities
rOpt = 0; % random cell r0
sdOpt = 0; % random cell self discharge 
cOpt = 0; % random cell-coulombic-efficiency values
lOpt = 0; % random cell leakage currents
rOpts = [tOpts, qOpt, rOpt, sdOpt, cOpt, lOpt];


packData = simRandPack(Ns,Nc,cycleFiles{cycleInd},model,rOpts);


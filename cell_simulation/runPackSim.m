% compile: mcc -I helper_function\ -I data\ -m runPackSim.m

% pack = runPackSim("3","3",'/data/drive_cycle_profiles/P14_us06_profile.mat',"/data/cell_models/P14model.mat", "1", "passive", "all",'[0,1,1,1,1,1]',"30", "[2, 2]") 

% pack = runPackSim("3","3",'data/drive_cycle_profiles/P14_us06_profile.mat',"data/cell_models/P14model.mat", "1", "active", "all",'[0,1,1,1,1,1]',"30", "[2, 2]")

function packData = runPackSim(Ns, Nc, cycleFile, cellModel, seed, balancing, sendSOCsWhen, randOps, sampleFactor, usageArray)

    filename = "simCell";

    % Run random pack simulaitons
    
    % Ns Number of cells
    % Nc Number of cycles to simulate
    Ns = str2double(Ns);
    Nc = str2double(Nc);

    % Downsampling, will sample every sampleFactor samplesC:\Users\A493192\Downloads\EVsim
    % When sampleFactor is 0, all data points will be sampled and sent
    sampleFactor = str2double(sampleFactor);

    randOps = strrep(randOps,","," ");
    randOps = str2num(randOps); %#ok<ST2NM> 
    
    % Percentage of time per day that the battery is under use (charging
    % and discharging). 1 - utilizationPercent is percentage of time per
    % day the battery is under resting conditions.
    usageArray = strrep(usageArray,","," ");
    usageArray = str2num(usageArray); %#ok<ST2NM> 

    seed = str2double(seed);
 
    if ~isdeployed; addpath data; end
    
    % "..\data\P14model.mat"
    load(cellModel, "model"); % Data Obtained from Dynamic Processing
    cycleProfile = load(cycleFile);
   
    % Dynamic model with milti temperature feature and different cell models
    % fname = 'A123modeldyn.json'; 
    % fid = fopen(fname); 
    % raw = fread(fid,inf); 
    % str = char(raw'); 
    % fclose(fid); 
    % model = jsondecode(str);
    
    % Enable/Disble randomized variables
%     tOpts = 0; % random cell teperatures
%     qOpt = 1; % random cell capacities
%     rOpt = 1; % random cell r0
%     sdOpt = 1; % random cell self discharge 
%     cOpt = 1; % random cell-coulombic-efficiency values
%     lOpt = 1; % random cell leakage currents
%     randOps = [tOpts, qOpt, rOpt, sdOpt, cOpt, lOpt];
    
    packData = simRandPack(Ns,Nc,cycleProfile,model,seed,balancing,sendSOCsWhen,randOps,filename,sampleFactor,usageArray);

%     f = figure();
%     plot(packData.storez');
%     waitfor(f);

end


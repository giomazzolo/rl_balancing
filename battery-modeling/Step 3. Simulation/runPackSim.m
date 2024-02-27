% compile: mcc -I ..\helper_function\ -I ..\data\ -m runPackSim.m

% pack = runPackSim("3","30",'..\data\uddsPower.mat',"..\data\P14model_dynamic.mat",'[0,1,1,1,1,1]')

function packData = runPackSim(Ns,Nc,cycleFile,cellModel,randOps, sample_rate)

    filename = "simCell";

    % Run random pack simulaitons
    
    % Ns Number of cells
    % Nc Number of cycles to simulate
    Ns = str2double(Ns);
    Nc = str2double(Nc);

    % Downsampling, will sample every sample_rate cycles
    % When sample_rate is 0, all data points will be sampled and sent
    sample_rate = str2double(sample_rate);

    randOps = strrep(randOps,","," ");
    randOps = str2num(randOps); %#ok<ST2NM> 

    % cycleFile = {'nyccPower.mat','uddsPower.mat','us06Power.mat','hwfetPower.mat'}; % drive cycles power demand files
 
    if ~isdeployed; addpath ..\data; end
    
    % "..\data\P14model_dynamic.mat"
    load(cellModel, "model"); % Data Obtained from Dynamic Processing
    
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
    
    packData = simRandPack(Ns,Nc,cycleFile,model,randOps,filename,sample_rate);

%     f = figure();
%     plot(packData.storez');
%     waitfor(f);

end


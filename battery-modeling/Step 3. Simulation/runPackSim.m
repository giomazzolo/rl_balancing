% compile: mcc -I ..\helper_function\ -I ..\data\ -m runPackSim.m

% pack = runPackSim("3","30",'uddsPower.mat',"..\data\P14model_dynamic.mat",'[0,1,1,1,1,1]')

function packData = runPackSim(Ns,Nc,cycleFile,cellModel,randOps)

    % Run random pack simulaitons
    
    % Ns Number of cells
    % Nc Number of cycles to simulate
    Ns = str2double(Ns);
    Nc = str2double(Nc);

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
    
    packData = simRandPack(Ns,Nc,cycleFile,model,randOps);

%     f = figure();
%     plot(packData.storez');
%     waitfor(f);

end


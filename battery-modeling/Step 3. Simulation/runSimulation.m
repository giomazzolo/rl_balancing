% Compile: mcc -I ..\helper_function\ -m runSimulation.m

function runSimulation()
    
    send2Py = false;
    plotSim_f = true;
    % close;clear;clc;
    
    if ~isdeployed
        addpath ..\helper_function\
    end
    
    load ..\data\P14_DYN_50_P25.mat DYNData
    load ..\data\P14model_dynamic.mat model % Data Obtained from Dynamic Processing
    
    deltaT = 1; 
    time = DYNData.script1.time - DYNData.script1.time(1);    
    t = (0:deltaT:time(end));
    voltage = interp1(time,DYNData.script1.voltage,t);
    current = interp1(time,DYNData.script1.current,t);
    time = t;
    
    if send2Py == true
        % Simulate and send
        %[vest,rck,hk,zk,sik,OCV] = simCell2py_tcp(current,35,deltaT,model,1,0,0, 40404);
        %[vest,rck,hk,zk,sik,OCV] = simCell2py_udp(current,35,deltaT,model,1,0,0, 40404, 40408);
        [vest,rck,hk,zk,sik,OCV] = simCell2py_mmap(current,25,deltaT,model,1,0,0);
    else
        % Simulate locally
        [vest,rck,hk,zk,sik,OCV] = simCellStep(current,25,deltaT,model,1,0,0);
    end
    
    if plotSim_f == true
        % for visualization purposes, plot the measured and simulated voltage data
        %subplot(1,2,1)
        %plot(time/3600,vest); % factor of 3600 converts seconds -> hours
        plot(time/3600,vest,time/3600,voltage); % factor of 3600 converts seconds -> hours
        %xlabel('Time (hr)'); ylabel('Voltage (V)'); title('Comparing measured to simulated voltage');
        %legend('Measured voltage','Simulated voltage');
        
        % Now, plot the voltage prediction error
        %subplot(1,2,2)
        %plot(time/3600,1000*(voltage-vest'));
        %xlabel('Time (hr)'); ylabel('Voltage (mV)'); title('Voltage prediction error');
    end
   
end
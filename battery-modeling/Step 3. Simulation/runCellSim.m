close;clear;clc;

send2Py = true;
plotData = false;

load ..\data\P14_DYN_50_P35.mat
load ..\data\P14model_dynamic.mat % Data Obtained from Dynamic Processing

deltaT = 1; 
time = DYNData.script1.time - DYNData.script1.time(1);    
t = (0:deltaT:time(end));
voltage = interp1(time,DYNData.script1.voltage,t);
current = interp1(time,DYNData.script1.current,t);
time = t;

if send2Py == true
    [vest,rck,hk,zk,sik,OCV] = simCell2py_mmap(current,35,deltaT,model,1,0,0);
else
    % Simulate locally
    [vest,rck,hk,zk,sik,OCV] = simCellStep(current,35,deltaT,model,1,0,0);
end

if plotData == true
    % for visualization purposes, plot the measured and simulated voltage data
    subplot(1,2,1)
    plot(time/3600,vest,time/3600,voltage); % factor of 3600 converts seconds -> hours
    xlabel('Time (hr)'); ylabel('Voltage (V)'); title('Comparing measured to simulated voltage');
    legend('Measured voltage','Simulated voltage');
    
    % Now, plot the voltage prediction error
    subplot(1,2,2)
    plot(time/3600,1000*(voltage-vest'));
    xlabel('Time (hr)'); ylabel('Voltage (mV)'); title('Voltage prediction error');
end
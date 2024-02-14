% function [vk,rck,hk,zk,sik,OCV] = simCell(ik,T,deltaT,model,z0,iR0,h0)
% ik - current, where (+) is discharge
% T  - temperature (degC)
% deltaT = sampling interval in data (s)
% model - standard model structure
% z0 - initial SOC
% iR0 - initial resistor currents as column vector
% h0 - initial hysteresis state

function [vk,rck,hk,zk,sik,OCV] = simCell2py_mmap(ik,T,deltaT,model,z0,iR0,h0)

    m_in = memmapfile('simCell_mmap_out.dat','Format','double');
   
    while m_in.Data(1) == 0; end % blocking wait

    m_out = memmapfile('simCell_mmap_in.dat','Format','double');
    m_out.Writable = true;
    
    disp("Starting...")
    
    % Define structure to send to python as json
   
    send2Py.state = true;
    send2Py.vk = 0;
    send2Py.sync = m_in.Data(1);
    
    % Force data to be column vector(s)
    % addpath ..\helper_function\
    ik = ik(:); iR0 = iR0(:);
    
    % Get model parameters from model structure
    RCfact = exp(-deltaT./abs(getParamESC('RCParam',T,model)))';
    G = getParamESC('GParam',T,model);
    Q = getParamESC('QParam',T,model);
    M = getParamESC('MParam',T,model);
    M0 = getParamESC('M0Param',T,model);
    RParam = getParamESC('RParam',T,model);
    R0Param = getParamESC('R0Param',T,model);
    etaParam = getParamESC('etaParam',T,model); % Coulombic Efficency
    
    etaik = zeros(1,length(ik));
    if ik(1) < 0
        etaik(1) = etaParam*ik(1); % Multiplying charging current with Coulombic Efficency
    end
    
    %%% First Iteration

    % Simulate the dynamic states of the model
    rck = zeros(length(RCfact),length(etaik));
    rck(:,1) = iR0; %Current through parallel resistor.
    
    csum = 0;
    zk = zeros(1,length(ik));
    zk(1) = z0-csum*deltaT/(Q*3600); % SOC
    
    if any(zk(1)>1.1)
        warning('Current may have wrong sign as SOC > 110%');
    end
    
    % Hysteresis stuff
    hk=zeros([length(ik) 1]);
    fac = hk*0;

    hk(1) = h0; sik = 0*hk;
    fac(1) = exp(-abs(G*etaik(1)*deltaT/(3600*Q)));
    sik(1) = sign(ik(1));
    
    % Compute output equation
    OCV = OCVfromSOCtemp(zk(1),T,model);

    vk = zeros(length(ik),1);
    vk(1) = OCV - rck(1)*RParam - ik(1)*R0Param + M*hk(1) + M0*sik(1); %Voltage Estimation
    send2Py.zk = zk(1);
    
    m_out.Data(3) = send2Py.zk;
    m_out.Data(2) = send2Py.state;
    m_out.Data(1) = send2Py.sync;

    %%% Remaining iterations
    for s = 2:length(ik)
        
        while m_in.Data(1) == send2Py.sync; end % blocking wait

        send2Py.sync = m_in.Data(1);

        if s == length(ik)-1
            send2Py.state = false;
        end
   
        if ik(s) < 0
            etaik(s) = etaParam*ik(s); % Multiplying charging current with Coulombic Efficency
        else
            etaik(s) = ik(s);
        end
        
        % Simulate the dynamic states of the model
        rck(s) = diag(RCfact)*rck(s-1) + (1-RCfact)*etaik(s-1);
        
        csum = csum + etaik(s-1);
        zk(s) = z0-csum*deltaT/(Q*3600); % SOC
        
        if zk(s)>1.1
            warning('Current may have wrong sign as SOC > 110%');
        end
        
        % Hysteresis stuff
        fac(s) = exp(-abs(G*etaik(s)*deltaT/(3600*Q)));
        hk(s)=fac(s-1)*hk(s-1)+(fac(s-1)-1)*sign(ik(s-1));
        sik(s) = sign(ik(s));
        if abs(ik(s))<Q/100, sik(s) = sik(s-1); end
        
        % Compute output equation
        OCV = OCVfromSOCtemp(zk(s),T,model);
        
        vk(s) = OCV - rck(s)*RParam - ik(s)*R0Param + M*hk(s) + M0*sik(s); %Voltage Estimation
        send2Py.zk = zk(s);

        m_out.Data(3) = send2Py.zk;
        m_out.Data(2) = send2Py.state;
        m_out.Data(1) = send2Py.sync;
        
    end

    clear m_out;
    clear m_in;
   
return
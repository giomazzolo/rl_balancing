% --------------------------------------------------------------------
% simRandPack: Simulate battery pack having Ns cells in series for Nc
% discharge/charge cycles, where all cells in pack can have random
% parameter values (e.g., capacity, resistance, etc.)
%
% Assumes no hysteresis in the cell model (this could be changed
% fairly easily; hysteresis makes results more difficult to interpret,
% so this assumption is okay for a first analysis, at least).
%
% Comunicates with python and sends data via the memory map data
% interface.
% --------------------------------------------------------------------


function packData = simRandPack(Ns, Nc, cycleFile, model, randOptions, filename, sampleFactor, utilizationPercent)

% State machine variable to configure the simulation by steps and allow
% resetting of the simulation in real time
simState = 0; % Initial state

% End simulation flag
endSim = false;

while endSim ~= true

switch simState
    
    case 0
    %% Initialize memmap
        simState = 1; % Next state

        % Memory map init
        m_in = memmapfile(strcat(filename,'_mmap_out.dat'),'Format','double');
        m_out = memmapfile(strcat(filename,'_mmap_in.dat'),'Format','double');
        m_out.Writable = true;
        send2Py.sync = 1.0;
        send2Py.state = true;
        outsize = length(m_out.Data)-2;

    case 1
    %% One time-initialized variables
        simState = 2; % Next state

        tOpt = randOptions(1); qOpt = randOptions(2); rOpt = randOptions(3);
        sdOpt = randOptions(4); cOpt = randOptions(5); lOpt = randOptions(6);
        profile = load(cycleFile, '-ascii'); % e.g., 'uddsPower.txt'
        
        if ~isdeployed; addpath ..\helper_function\; end
        
        % ------------------------------------------------------------------
        % Create storage for all cell states after completion of each cycle
        % ------------------------------------------------------------------
        if ~isdeployed  % Dont store anything if compiled
            packData.storez = zeros([Ns Nc]); % create storage for final SOC
            packData.storeirc = zeros([Ns Nc]);
        end
        % ------------------------------------------------------------------
        % Initialize default states for ESC cell model
        % ------------------------------------------------------------------
        maxSOC = 0.85; % cell SOC when pack is "fully charged"
        minSOC = 0.15; % cell SOC when pack is "fully discharged"

        utilizationInSec = utilizationPercent * 24 * 60 * 60;
        restingInSec = (1-utilizationPercent) * 24 * 60 * 60;
        
    case 2
    %% Initialize variables
        simState = 3; % Next state

        % Downsampling factor counter
        sample_cnt = 0;        
        %Balancing parameters
        rl_balance_i = zeros([Ns 1]);

        z = maxSOC*ones(Ns,1); % start fully charged
        irc = zeros(Ns,1); % at rest
        ik = zeros([Ns 1]); % current experienced by each cell
        
        % TODO: Implement multi-temperature model compatibility %
        % WORKAROUND: Force temperature to be the same for all simulated cells
        % until multi-cell temperature is implemented. This might be
        % skipped and left as-is for simplicity.
        % Set cell temperatures based on tOpt
        if tOpt % set to "if 1," to execute, or "if 0," to skip this code
            %T = 22.5 + 5*rand([Ns 1]);
            T = 22.5 + 5*rand(1);
        else
            %T = 25*ones([Ns 1]);
            T = 25;
        end
        
        % Set self-discharge "cell temperature"
        Tsd = T - 5 + 10*rand([Ns 1]);
        % Set cell module leakage current based on lOpt
        
        if lOpt
            leak = 0.01 + 0.002*rand([Ns 1]);
        else
            leak = 0.01*ones([Ns 1]);
        end
        % ------------------------------------------------------------------
        % Default initialization for cells within the pack
        % Note that since T has Ns elements, there is one parameter value
        % per cell (even if all turn out to be identical)
        % ------------------------------------------------------------------
        q = getParamESC('QParam',T,model);
        rc = exp(-1./abs(getParamESC('RCParam',T,model)));
        r = (getParamESC('RParam',T,model)).*(1-rc);
        r0 = getParamESC('R0Param',T,model);
        rt = 2*0.000125; % 125 microOhm resistance for each tab
        maxVlim = OCVfromSOCtemp(maxSOC,T,model);
        minVlim = OCVfromSOCtemp(minSOC,T,model);
        eta = ones([Ns 1]);
        % ------------------------------------------------------------------
        % Modified initialization for cell variability
        % ------------------------------------------------------------------
        % Set individual random cell-capacity values
        if qOpt % set to "if 1," to execute, or "if 0," to skip this code
            q=q-0.25+0.5*rand([Ns 1]); % random capacity for ea. cell
        end
        % Set individual random cell-resistance values
        if rOpt % set to "if 1," to execute, or "if 0," to skip this code
            r0 = r0-0.0005+0.0015*rand(Ns,1);
        end
        r0 = r0 + rt; % add tab resistance to cell resistance
        R = sum(r0,1);
        % Set individual random cell-coulombic-efficiency values
        if cOpt % set to "if 1," to execute, or "if 0," to skip this code
            eta = eta - 0.001 - 0.002*rand([Ns 1]);
        end

    case 3 
    %% Cell pack simulation
        % Now, simulate pack performance using ESC cell model.
        % ------------------------------------------------------------------
        theCycle = 1; theState = 'discharge';
        disCnt = 0; % start at beginning of profile
        fprintd(' Cycle = 1, discharging... ');

        simState = 4; % Next state

        utilizationCounter = 0;
        restingCounter = 0;
        
        while theCycle <= Nc
            
            v = OCVfromSOCtemp(z,T,model); % get OCV for each cell
            v = v - r.*irc; % add in capacitor voltages
            V = sum(v); % Total voltage excluding I*R
            vt = v-ik.*r0; % Cell terminal voltages
            % Hystheresis can be added here

            switch( theState )
                case 'discharge'
                    % Get instantaneous demanded pack power, repeating profile
                    P = profile(rem(disCnt,length(profile))+1);
                    % Compute demanded pack current based on unloaded voltage
                    I = V/(2*R) - sqrt(V^2/R^2 - 4*P/R)/2;
                    % Default cell current = pack current
                    ik = I*ones(Ns,1);
                    if I < 0 % If we happen to be charging this momement
                        ik = ik.*eta;
                    end
                    if min(z) <= minSOC || min(vt) < minVlim % stop discharging
                        theState = 'charge';
                        chargeFactor = 1;
                        ik = 0*ik;
                        fprintd('charging... ');
                    end

                    if utilizationCounter >= utilizationInSec
                        fprintd('resting... ');
                        theState = 'resting';
                        utilizationCounter = 0;
                    end

                    disCnt = disCnt + 1;
                    utilizationCounter = utilizationCounter + 1;

                case 'charge'
                    % start charging @ 6.6kW, then taper
                    P = -6600/chargeFactor;
                    I = V/(2*R) - sqrt(V^2/R^2 - 4*P/R)/2;
                    I = max(-min(q),I); % limit to 1C charge rate max
                    ik = I*eta; % Charge coulombic eff.
                    if max(vt)>=maxVlim
                        if chargeFactor > 32 % bail after 6.6kW/32 charge
        
                            if ~isdeployed % Dont store anything if compiled
                                packData.storez(:,theCycle) = z;
                                packData.storeirc(:,theCycle) = irc;
                            end
        
                            theState = 'discharge';
                            disCnt = 0;
                            ik = 0*ik;
                            theCycle = theCycle + 1;
                            if theCycle <= Nc
                                fprintd('\n Cycle = %d, discharging... ',theCycle);
                            end
                        end
                        chargeFactor = chargeFactor*2;
                    end

                case 'resting'

                    % Balancing from RL model feedback
                    ik = rl_balance_i;

                    if restingCounter >= restingInSec
                        fprintd('discharging... ');
                        theState = 'discharge';
                        restingCounter = 0;
                    end

                    restingCounter = restingCounter + 1;
                    sample_cnt = sample_cnt + 1;

                otherwise
                    error('charge/discharge state has been corrupted')
            end
            % Simulate self discharge via variable resistor in parallel
            if sdOpt == 1
                rsd = ((-20+0.4*Tsd).*z + (35-0.5*Tsd))*1e3;
                ik = ik + vt./rsd;
            end
            % Simulate leakage current
            ik = ik + leak;
                    
            % Update each cell SOC
            z = z - (1/3600)*ik./q;
            % Update resistor currents
            irc = rc.*irc + (1-rc).*ik;

            %% Balancing Data sent and recieved here
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if sample_cnt >= sampleFactor % Downsampling
                sample_cnt = 0;
        
                % Write data to memmap file        
                m_out.Data(3:end) = z(1:outsize);
                m_out.Data(2) = send2Py.state;
                m_out.Data(1) = send2Py.sync;
        
                % Memory map sync
                while m_in.Data(1) ~= send2Py.sync; end % blocking wait
                send2Py.sync = m_in.Data(1) + 1;
                cmd = m_in.Data(2);
        
                if cmd == 1 % Stop simulation prematurely
                    simState = 4;
                    break;        
                elseif cmd == 2 % Reset simulation prematurely
                    % Go back to variable init state and re-initialize the
                    % simulation
                    simState = 2;
                    break;
                end
        
                % Get balancing feedback from the RL mode
                % Feedback is trearted as a discharge current value
                rl_balance_i = m_in.Data(3:end);
         
            end
        
            %sample_cnt = sample_cnt + 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        end % end while

    case 4
    % Finalize simulation
       
        % Communicate that somulation has finished
        send2Py.state = false;
        %m_out.Data(3:end) = z(1:outsize);
        m_out.Data(2) = send2Py.state;
        m_out.Data(1) = send2Py.sync;

        % Wait for command
        while m_in.Data(1) ~= send2Py.sync; end % blocking wait
        send2Py.sync = m_in.Data(1) + 1;
        cmd = m_in.Data(2);

        if cmd == 1 % Stop simulation

            endSim = true;
            clear m_out;
            clear m_in;
            
            fprintd('\n');
            
            if ~isdeployed  % Dont store anything if compiled
                packData.q = q; packData.rc = rc; packData.eta = eta;
                packData.r = r; packData.r0 = r0; packData.Tsd = Tsd;
                packData.T = T; packData.leak = leak;
            end

        elseif cmd == 2 % Reset simulation
            % Go back to variable init state and re-initialize the
            % simulation
            send2Py.state = true;
            simState = 2;
        else
            simState = -1;
        end

    otherwise % Should be unreachable
        endSim = true;
end
end
end

%% Local functions

% Debug print function
% Only prints if run in script form
% Does not print if script is compiled
function fprintd(varargin)
    if ~isdeployed
        if nargin==2
            fprintf(varargin{1},varargin{1});
        else
            fprintf(varargin{1});
        end
    end
end


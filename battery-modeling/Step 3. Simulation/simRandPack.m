% --------------------------------------------------------------------
% simRandPack: Simulate battery pack having Ns cells in series for Nc
% discharge/charge cycles, where all cells in pack can have random
% parameter values (e.g., capacity, resistance, etc.)
%
% Assumes no hysteresis in the cell model (this could be changed
% fairly easily; hysteresis makes results more difficult to interpret,
% so this assumption is okay for a first analysis, at least).
% --------------------------------------------------------------------

function packData = simRandPack(Ns,Nc,cycleFile,model,randOptions)

tOpt = randOptions(1); qOpt = randOptions(2); rOpt = randOptions(3);
sdOpt = randOptions(4); cOpt = randOptions(5); lOpt = randOptions(6);
profile = load(cycleFile); % e.g., 'uddsPower.txt'

% ------------------------------------------------------------------
% Create storage for all cell states after completion of each cycle
% ------------------------------------------------------------------
packData.storez = zeros([Ns Nc]); % create storage for final SOC
packData.storeirc = zeros([Ns Nc]);
% ------------------------------------------------------------------
% Initialize default states for ESC cell model
% ------------------------------------------------------------------
maxSOC = 0.95; % cell SOC when pack is "fully charged"
minSOC = 0.1; % cell SOC when pack is "fully discharged"
z = maxSOC*ones(Ns,1); % start fully charged
irc = zeros(Ns,1); % at rest
ik = zeros([Ns 1]); % current experienced by each cell

% Set cell temperatures based on tOpt
if tOpt % set to "if 1," to execute, or "if 0," to skip this code
    T = 22.5 + 5*rand([Ns 1]);
else
    T = 25*ones([Ns 1]);
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

% Now, simulate pack performance using ESC cell model.
% ------------------------------------------------------------------
theCycle = 1; theState = 'discharge';
disCnt = 0; % start at beginning of profile
fprintf(' Cycle = 1, discharging... ');
while theCycle <= Nc
    v = OCVfromSOCtemp(z,T,model); % get OCV for each cell
    v = v - r.*irc; % add in capacitor voltages
    V = sum(v); % Total voltage excluding I*R
    vt = v-ik.*r0; % Cell terminal voltages
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
            if min(z) <= minSOC || min(vt) < minVlim, % stop discharging
                theState = 'charge';
                chargeFactor = 1;
                ik = 0*ik;
                fprintf('charging... ');
            end
            disCnt = disCnt + 1;
        case 'charge'
            % start charging @ 6.6kW, then taper
            P = -6600/chargeFactor;
            I = V/(2*R) - sqrt(V^2/R^2 - 4*P/R)/2;
            I = max(-min(q),I); % limit to 1C charge rate max
            ik = I*eta; % Charge coulombic eff.
            if max(vt)>=maxVlim,
                if chargeFactor > 32, % bail after 6.6kW/32 charge
                    packData.storez(:,theCycle) = z;
                    packData.storeirc(:,theCycle) = irc;
                    theState = 'discharge';
                    disCnt = 0;
                    ik = 0*ik;
                    theCycle = theCycle + 1;
                    if theCycle <= Nc,vbvggf
                        fprintf('\n Cycle = %d, discharging... ',theCycle);
                    end
                end
                chargeFactor = chargeFactor*2;
            end
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
    z = z - (1/3600)*ik./q; % Update each cell SOC
    irc = rc.*irc + (1-rc).*ik; % Update resistor currents
end % end while

fprintf('\n');
packData.q = q; packData.rc = rc; packData.eta = eta;
packData.r = r; packData.r0 = r0; packData.Tsd = Tsd;
packData.T = T; packData.leak = leak;

end

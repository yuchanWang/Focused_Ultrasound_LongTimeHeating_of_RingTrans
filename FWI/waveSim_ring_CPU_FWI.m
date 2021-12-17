function [simData, gradData, tArr, rcvGridIdx] = waveSim_ring_CPU_FWI(numElements, medium, cellSize, imgSize, adjoint, normalize, waveDataExp,sensorX,sensorY,sensor_freq)

%% ========================================================================
% DEFINITIONS
%
% Simulates wave fields given acoustical image and transducer locations
% 
% INPUTS
% numElements -> number of transducer elements in total
% medium -> medium used in k-Wave (SOS + attenuation images)
% cellSize -> physical grid spacing of the image
% imgSize -> image size in pixels
% adjoint -> 0: wave fields simulation, 1: adjoint fields simulation (dirichlet boundary condition enforced)
% normalize -> 0: do not normalize the gradient, 1: normalize the gradient
% waveDiff -> input rasidual wave fields used in adjoint mode
% 
% OUTPUTS
% simData -> simulated wave fields at all pixel locations
% gradData -> gradient vector
% tArr -> time array
% rcvGridIdx -> pixel locations of all receivers
%
% AUTHORS
% Rungroj Jintamethasawat (rungroj@umich.edu) and Yunhao Zhu (zyh6557@126.com)
% =========================================================================

%% ========================================================================
% INITIALIZATION
% =========================================================================

% adjoint = 0;
% normalize = 0;
% waveDataExp = 0;

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 10;
PML_Y_SIZE = 10;

% create the computational grid
Nx = imgSize(1);        % number of grid points in the x direction
Ny = imgSize(2);        % number of grid points in the y direction
dx = cellSize;
dy = cellSize;

% create the simulation grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% create the time array
dt  = 1 / sensor_freq;  
t_end = sqrt( kgrid.x_size.^2 + kgrid.y_size.^2 ) / 1450;
Nt = round(t_end / dt);
% create the time array
kgrid.setTime(Nt, dt);

C = zeros(Nx, Ny);
for i = 1 : numElements
    C(sensorX(i), sensorY(i)) = 1;
end
% get the pixels corresponding to receivers
sensorVolume = C;
rcvGridIdx = find(sensorVolume == 1);
sensor.mask = ones(Nx, Ny);

% source for each element
waveSource = cell(1, numElements);

% define the properties of the pulse used to drive the transducer
tone_burst_freq = 0.5E6; % (Hz)
tc = t_end/8;
devia = 0.8e-6;

% define input signal pressure (see K. Wang's paper)
input_signal = zeros(1, length(kgrid.Nt));

for i = 1:kgrid.Nt
    t = kgrid.dt*i;
    input_signal(i) = exp(-(t-tc)^2/(2*devia^2))*sin(2*pi*tone_burst_freq*t);
end

for i = 1:numElements
    waveSource{1, i} = input_signal;
end

%% ========================================================================
% WAVE FIELD SIMULATION
% =========================================================================

% number of simulated time points 
numTimePoints = length(kgrid.t_array);

% simulation data
simData = zeros(numTimePoints, numElements, numElements, 'single');

% parameters for simulation
parameterCell = cell(1, numElements);

% for each element, compute the simulated waveform
for acqNum = 1:numElements
    source_sim.p_mask = zeros(Nx, Ny);

    source_sim.p_mask(rcvGridIdx(acqNum)) = 1; % activate only one element
    source_sim.p = waveSource{1, acqNum}; % initialize the source pressure with predefined pressure
    
    % assign the input options (see k-Wave documentation for more details)
    % Note: change PlotSim to 'true' if wanting to show the simulation
    input_args = {'DisplayMask', source_sim.p_mask, 'PMLInside', false, 'PlotPML', false, ...
    'PMLSize', [PML_X_SIZE, PML_Y_SIZE], 'PlotSim', false, 'RecordMovie', false, 'Smooth', [true, true, true], ...
    'DataCast', 'single'};

    % pack the necessary data into the parameter cell
    parameterCell{1, acqNum} = {kgrid, medium, source_sim, sensor, input_args, acqNum};
end

%% ========================================================================
% SIMULATION
% =========================================================================

gradData = zeros(numElements, Nx*Ny);

gradDataNorm = zeros(numElements, Nx*Ny);

fprintf('Running k-wave simulation...\n');

% simulate the transducer field patterns

if(~adjoint)
    waveDataExp = simData;
end

parfor k = 1:numElements
    parameterList = parameterCell{1, k};
    kgridParam = parameterList{1, 1};
    mediumParam = parameterList{1, 2};
    sourceParam = parameterList{1, 3};
    sensorParam = parameterList{1, 4};
    inputArgsParam = parameterList{1, 5};
    elNumParam = parameterList{1, 6};

    [simData(:, :, k), gradData(k, :), gradDataNorm(k, :)] = calculationSlave_CPU_FWI(kgridParam, mediumParam, [PML_X_SIZE, PML_Y_SIZE], ...
        sourceParam, sensorParam, inputArgsParam, rcvGridIdx, adjoint, waveDataExp(:, :, k), kgrid.dt, elNumParam);

end

if(adjoint)
    gradData = sum(gradData, 1);
    
    % cancel out the geometric spreading
    if(normalize)
        gamma = 10^-12;
        normJ = zeros(1, Nx*Ny);

        for m = 1:numElements
            normJ = normJ + gradDataNorm(m, :);
        end
        normJ = sqrt(normJ + gamma^2);

        gradData = gradData./normJ;
    end
end

tArr = kgrid.t_array;
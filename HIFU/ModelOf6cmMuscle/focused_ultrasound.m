% % Heating By A Focused Ultrasound Transducer
% %
% % This example demonstrates how to combine acoustic and thermal simulations
% % in k-Wave to calculate the heating by a focused ultrasound transducer. It
% % builds on the Simulating Transducer Field Patterns and Using A Binary
% % Sensor Mask examples.
% % 
% % author: Bradley Treeby
% % date: 27th April 2017
% % last update: 28th April 2017
% %  
% % This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% % Copyright (C) 2017 Bradley Treeby
% 
% % This file is part of k-Wave. k-Wave is free software: you can
% % redistribute it and/or modify it under the terms of the GNU Lesser
% % General Public License as published by the Free Software Foundation,
% % either version 3 of the License, or (at your option) any later version.
% % 
% % k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% % WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% % FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% % more details. 
% % 
% % You should have received a copy of the GNU Lesser General Public License
% % along with k-Wave. If not, see <http://www.gnu.org/licenses/>.
% 
clearvars;
% =========================================================================
% ACOUSTIC SIMULATION
% =========================================================================

% define the PML size
pml_size = 10;              % [grid points]

% define the grid parameters
Nx = 510 - 2 * pml_size;    % [grid points]
Ny = 510 - 2 * pml_size;    % [grid points]
Nz = 60 - 2 * pml_size;     %300/988
dx = 220e-3/Nx;               % [m]
dy = 220e-3/Ny;               % [m]
dz = dx;                      % [m]
% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
muscleRadius = 6e-2;
fatRadius = muscleRadius + 1e-2;
skinRadius = fatRadius + 1.7e-3;
transRadius = 10e-2;
tumorRadius = 6.5e-3;
boneRadius = 1.25e-2;
boneThickness = 0.58e-2;

muscleSpeed = 1580;
waterSpeed = 1482;%1482.3  % 4.7061e+05   1510;   4.7327e+05   1482;
fatSpeed = 1430;
marrowSpeed = 1371.9;
skinSpeed = 1595;
tumorSpeed = 1450;
boneSpeed = 2198;
boneDistance = 30e-3;
tumorDistance = 25e-3;

muscleAtten = 0.57;
fatAtten = 0.6;
waterAtten = 2.17e-3;
skinAtten = 0.264;
boneAtten = 2.54;
tumorAtten = 0.4; %[dB/(MHz^y cm)]

muscleDensity = 1090;
waterDensity  = 1000;
fatDensity    = 911;
skinDensity   = 1109;
% tumorDensity  = 1079;
tumorDensity  = 1050;
boneDensity   = 1178;

bloodSpeed = 1550;
bloodAtten = 0.3;
bloodDensity = 1050;

vesselWallSpeed = 1570;
vesselWallAtten = 0.54;
vesselWallDensity = 1102;

% define medium properties related to perfusion
vesselRadius1    = 3.4e-3;
vesselRadius2    = 1.15e-3;
vesselRadius3    = vesselRadius2;
vesselThickness  = 0.5e-3;

vesselDistance1    = 1e-2;
vesselDistance2    = 3.8e-2;
vesselDistance3    = 4.3e-2;
vesselDistancex4   = 3.2e-2;
%--------------------------------------------------------------%
% define the Density of the propagation medium
%--------------------------------------------------------------%
medium.sound_speed = waterSpeed*ones(Nx, Ny);  % [m/s]

a = makeDisc(Nx, Ny, Nx/2, Ny/2, muscleRadius/dx);
medium.sound_speed(a == 1) = muscleSpeed; 

b = 2*makeDisc(Nx, Ny, Nx/2, Ny/2, fatRadius/dx);
b(a == 1) = 1;
medium.sound_speed(b == 2) = fatSpeed;

c = 3*makeDisc(Nx, Ny, Nx/2, Ny/2, skinRadius/dx);
c(a==1) = 1;
c(b==2) = 2;
medium.sound_speed(c == 3) = skinSpeed;

f = 4*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, (boneRadius-boneThickness)/dx);
medium.sound_speed(f == 4) = marrowSpeed;
e = 5*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, boneRadius/dx);
e(f == 4) = 1;
medium.sound_speed(e == 5) = boneSpeed;

% d = makeBall(Nx, Ny, Nz, Nx/2+tumorDistance/dx, Ny/2, Nz/2, tumorRadius/dx);
d = makeDisc(Nx, Ny, Nx/2+(tumorDistance/dx), Ny/2, tumorRadius/dx);
medium.sound_speed(d == 1) = tumorSpeed;

g = 6*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance1/dx, vesselRadius1/dx);
ga = 7*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance1/dx, (vesselRadius1-vesselThickness)/dx)
g(ga == 7) = 1;
medium.sound_speed(ga == 7) = vesselWallSpeed;
medium.sound_speed(g == 6) = bloodSpeed;

h = 8*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance2/dx, vesselRadius2/dx); 
ha = 9*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance2/dx, (vesselRadius2-vesselThickness)/dx)
h(ha == 9) = 1;
medium.sound_speed(ha == 9) = vesselWallSpeed;
medium.sound_speed(h == 8) = bloodSpeed;

k = 10*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance3/dx, vesselRadius3/dx);  
ka = 11*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance3/dx, (vesselRadius3-vesselThickness)/dx)
k(ka == 11) = 1;
medium.sound_speed(ka == 11) = vesselWallSpeed;
medium.sound_speed(k == 10) = bloodSpeed;

medium.sound_speed = repmat(medium.sound_speed,[1,1,Nz]);
medium.sound_speed(Nx/2+(tumorDistance/dx)-20:Nx/2+(tumorDistance/dx)+20,Ny/2-20:Ny/2+20,1:Nz/2-7) = muscleSpeed;
medium.sound_speed(Nx/2+(tumorDistance/dx)-20:Nx/2+(tumorDistance/dx)+20,Ny/2-20:Ny/2+20,Nz/2+7:end) = muscleSpeed;

%--------------------------------------------------------------%
% define the density of the propagation medium
%--------------------------------------------------------------%
% medium.density = medium.sound_speed./1.5;
medium.density = waterDensity*ones(Nx, Ny);  % [m/s]

a = makeDisc(Nx, Ny, Nx/2, Ny/2, muscleRadius/dx);
medium.density(a == 1) = muscleDensity; 

b = 2*makeDisc(Nx, Ny, Nx/2, Ny/2, fatRadius/dx);
b(a == 1) = 1;
medium.density(b == 2) = fatDensity;

c = 3*makeDisc(Nx, Ny, Nx/2, Ny/2, skinRadius/dx);
c(a==1) = 1;
c(b==2) = 2;
medium.density(c == 3) = skinDensity;

f = 4*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, (boneRadius-boneThickness)/dx);
medium.density(f == 4) = fatDensity;
e = 5*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, boneRadius/dx);
e(f == 4) = 1;
medium.density(e == 5) = boneDensity;

d = makeDisc(Nx, Ny, Nx/2+(tumorDistance/dx), Ny/2, tumorRadius/dx);
medium.density(d == 1) = tumorDensity;

g = 6*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance1/dx, vesselRadius1/dx);
ga = 7*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance1/dx, (vesselRadius1-vesselThickness)/dx)
g(ga == 7) = 1;
medium.density(ga == 7) = vesselWallDensity;
medium.density(g == 6) = bloodDensity;

h = 8*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance2/dx, vesselRadius2/dx); 
ha = 9*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance2/dx, (vesselRadius2-vesselThickness)/dx)
h(ha == 9) = 1;
medium.density(ha == 9) = vesselWallDensity;
medium.density(h == 8) = bloodDensity;

k = 10*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance3/dx, vesselRadius3/dx);  
ka = 11*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance3/dx, (vesselRadius3-vesselThickness)/dx)
k(ka == 11) = 1;
medium.density(ka == 11) = vesselWallDensity;
medium.density(k == 10) = bloodDensity;

medium.density = repmat(medium.density,[1,1,Nz]);
medium.density(Nx/2+(tumorDistance/dx)-20:Nx/2+(tumorDistance/dx)+20,Ny/2-20:Ny/2+20,1:Nz/2-7) = muscleDensity;
medium.density(Nx/2+(tumorDistance/dx)-20:Nx/2+(tumorDistance/dx)+20,Ny/2-20:Ny/2+20,Nz/2+7:end) = muscleDensity;
%--------------------------------------------------------------%
% define attenuation 
%--------------------------------------------------------------%
medium.alpha_coeff = waterAtten*ones(Nx, Ny);  % [m/s]

a = makeDisc(Nx, Ny, Nx/2, Ny/2, muscleRadius/dx);
medium.alpha_coeff(a == 1) = muscleAtten; 

b = 2*makeDisc(Nx, Ny, Nx/2, Ny/2, fatRadius/dx);
b(a == 1) = 1;
medium.alpha_coeff(b == 2) = fatAtten;

c = 3*makeDisc(Nx, Ny, Nx/2, Ny/2, skinRadius/dx);
c(a==1) = 1;
c(b==2) = 2;
medium.alpha_coeff(c == 3) = skinAtten;

f = 4*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, (boneRadius-boneThickness)/dx);
medium.alpha_coeff(f == 4) = fatAtten;
e = 5*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, boneRadius/dx);
e(f == 4) = 1;
medium.alpha_coeff(e == 5) = boneAtten;

d = makeDisc(Nx, Ny, Nx/2+(tumorDistance/dx), Ny/2, tumorRadius/dx);
medium.alpha_coeff(d == 1) = tumorAtten;

g = 6*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance1/dx, vesselRadius1/dx);
ga = 7*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance1/dx, (vesselRadius1-vesselThickness)/dx)
g(ga == 7) = 1;
medium.alpha_coeff(ga == 7) = vesselWallAtten;
medium.alpha_coeff(g == 6) = bloodAtten;

h = 8*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance2/dx, vesselRadius2/dx); 
ha = 9*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance2/dx, (vesselRadius2-vesselThickness)/dx)
h(ha == 9) = 1;
medium.alpha_coeff(ha == 9) = vesselWallAtten;
medium.alpha_coeff(h == 8) = bloodAtten;

k = 10*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance3/dx, vesselRadius3/dx);  
ka = 11*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance3/dx, (vesselRadius3-vesselThickness)/dx)
k(ka == 11) = 1;
medium.alpha_coeff(ka == 11) = vesselWallAtten;
medium.alpha_coeff(k == 10) = bloodAtten;

medium.alpha_coeff = repmat(medium.alpha_coeff,[1,1,Nz]);
medium.alpha_coeff(Nx/2+(tumorDistance/dx)-20:Nx/2+(tumorDistance/dx)+20,Ny/2-20:Ny/2+20,1:Nz/2-7) = muscleAtten;
medium.alpha_coeff(Nx/2+(tumorDistance/dx)-20:Nx/2+(tumorDistance/dx)+20,Ny/2-20:Ny/2+20,Nz/2+7:end) = muscleAtten;

%--------------------------------------------------------------%
% define power 
%--------------------------------------------------------------%
medium.alpha_power = 1.5;

% calculate the time step using an integer number of points per period
freq = 1.5e6;             % [Hz]
amp = 50e3;           % [Pa] 150e3
ppw = muscleSpeed / (freq * dx); % points per wavelength
cfl = 0.3;                              % cfl number
ppp = ceil(ppw / cfl);                  % points per period
T   = 1 / freq;                         % period [s]
dt  = T / ppp;                          % time step [s]

% calculate the number of time steps to reach steady state
t_end = sqrt( kgrid.x_size.^2 + kgrid.y_size.^2 + kgrid.z_size.^2) / muscleSpeed; 
Nt = round(t_end / dt);

% create the time array
kgrid.setTime(Nt, dt);

% % define the input signal
% source.p = createCWSignals(kgrid.t_array, freq, amp, 0);


% set the sensor mask to cover the entire grid
sensor.mask = ones(Nx, Ny, Nz);
sensor.record = {'p'};

% record the last 3 cycles in steady state
num_periods = 10;
T_points = round(num_periods * T / kgrid.dt);
sensor.record_start_index = 1051 - T_points/2;
kgrid.Nt = 1051 + round(T_points/2);
% sensor.record_start_index = Nt - T_points + 1;

%% define the source parameters(set transducer)
source.p_mask = zeros(Nx,Ny,Nz);
load('Gridx_Pos.mat');
load('Gridy_Pos.mat');
num_trans = 256;
for i = Nz/2-10:2:Nz/2+10
    for j = 1:num_trans
        source.p_mask(x_Pos(j),y_Pos(j),i) = 1;
    end
end
elementX=[];
elementY=[];
elementZ=[];

for i = Nz/2-10:2:Nz/2+10
    [X,Y] = find(source.p_mask(:,:,i) == 1);
    Z = i*ones(size(X));
    elementX = [elementX;X];
    elementY = [elementY;Y];
    elementZ = [elementZ;Z];
end


%% calculate delay
sampling_freq = 1/kgrid.dt;     % [Hz]
tone_burst_freq = 1.5e6; 
tone_burst_cycles = 5;
furthest_distance = 278.2017*dx;
longest_time = furthest_distance/muscleSpeed;  %[s]
% p_total = zeros(Nx,Ny,Nz,T_points);
for i = 1281:1536
    tone_burst_offset(i) = longest_time-(dx*sqrt((elementX(i)-(Nx/2+tumorDistance/dx))^2+(elementY(i)-Ny/2)^2))/(muscleSpeed);
    for j = 1:5
        tone_burst_offset(i-256*j) = tone_burst_offset(i);
        tone_burst_offset(i+256*j) = tone_burst_offset(i);
    end
end
tone_burst_offset = tone_burst_offset/kgrid.dt;

% create the tone burst signals
source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, ...
    'SignalOffset', tone_burst_offset);
% num_ring_up_cycles = 0;
% num_ring_down_cycles = 0;
% source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles,'Envelope',[num_ring_up_cycles, num_ring_down_cycles],'SignalOffset', tone_burst_offset);

% source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles);
source_strength = 1.45e12;          % [Pa]
source.p = (source_strength ./ (muscleSpeed * muscleSpeed/1.5)) .* source.p;

num_source_time_points = length(source.p(1,:));

% get suitable scaling factor for plot axis
[~, scale, prefix] = scaleSI(kgrid.t_array(num_source_time_points));

% plot the input time series
% figure;
% stackedPlot(kgrid.t_array(1:num_source_time_points) * scale, source.p([1,round(end/2),end],:));
% xlabel(['Time [' prefix 's]']);
% ylabel('Input Signals');

% create a display mask to display the transducer
display_mask = source.p_mask;

% create a sensor mask covering the entire computational domain using the
% opposing corners of a rectangle
sensor.mask = [1, 1, 1, Nx, Ny, Nz].';

% set the record mode capture the final wave-field and the statistics at
% each sensor point 
% sensor.record = {'p_final', 'p_max', 'p_rms'};

% set the input arguements
input_args = {'PMLInside', false, 'PlotPML', false,'PlotSim', false, 'DisplayMask', ...
    'off', 'PlotScale', [-1, 1] * amp};

sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});

% =========================================================================
% CALCULATE HEATING
% =========================================================================

% convert the absorption coefficient to nepers/m
alpha_np = db2neper(medium.alpha_coeff, medium.alpha_power) * ...
    (2 * pi * freq).^medium.alpha_power;

% extract the pressure amplitude at each position
% p = extractAmpPhase(sensor_data.p, 1/kgrid.dt, freq);
p = extractAmpPhase(sensor_data.p, 1/kgrid.dt, freq);

% reshape the data, and calculate the volume rate of heat deposition
p = reshape(p, Nx, Ny, Nz);

Q = alpha_np .* p.^2 ./ (medium.density .* medium.sound_speed);
filename = 'AcousticFocusing.mat';save(filename);
% % =========================================================================
% % THERMAL SIMULATION
% % =========================================================================

%%
clearvars;
load('AcousticFocusing.mat');
clear medium source sensor;

% define medium properties related to diffusion
muscleDensity = 1090;
waterDensity  = 1000;
fatDensity    = 911;
skinDensity   = 1109;
% tumorDensity  = 1079;
tumorDensity  = 1050;
boneDensity   = 1178;

muscleTc = 0.49;
waterTc  = 0.6;  %0.6
fatTc    = 0.21;
skinTc   = 0.37;
tumorTc  = 0.52;
boneTc   = 0.32;

muscleSh = 3421;
waterSh  = 4184.4;  %4184.4
fatSh    = 2348;
skinSh   = 3391;
tumorSh  = 3540;
boneSh   = 1313;

% define medium properties related to perfusion
vesselRadius1    = 3.4e-3;
vesselRadius2    = 1.15e-3;
vesselRadius3    = vesselRadius2;


vesselDistance1    = 1e-2;
vesselDistance2    = 3.8e-2;
vesselDistance3    = 4.3e-2;
vesselDistancex4   = 3.2e-2;


bloodDensity                       = 1060;     %1060 [kg/m^3]
bloodSh                            = 3617;     % [J/(kg.K)]
bloodPr                            = 0.01;     % [1/s]
medium.blood_ambient_temperature   = 37;       % [degC]

%--------------------------------------------------------------%
% define the Density of the propagation medium
%--------------------------------------------------------------%
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

d = makeDisc(Nx, Ny, Nx/2+tumorDistance/dx, Ny/2, tumorRadius/dx);
medium.density(d == 1) = tumorDensity;

medium.density = repmat(medium.density,[1,1,Nz]);
medium.density(Nx/2+tumorDistance/dx-20:Nx/2+tumorDistance/dx+20,Ny/2-20:Ny/2+20,1:Nz/2-7) = muscleDensity;
medium.density(Nx/2+tumorDistance/dx-20:Nx/2+tumorDistance/dx+20,Ny/2-20:Ny/2+20,Nz/2+7:end) = muscleDensity;% A = medium.density;

medium.blood_density = zeros(Nx, Ny);
g = 6*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance1/dx, vesselRadius1/dx);  
medium.blood_density(g == 6) = bloodDensity-muscleDensity;  
h = 7*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance2/dx, vesselRadius2/dx);  
medium.blood_density(h == 7) = bloodDensity-muscleDensity;
k = 8*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance3/dx, vesselRadius3/dx);  
medium.blood_density(k == 8) = bloodDensity-muscleDensity;
medium.blood_density = repmat(medium.blood_density,[1,1,Nz]);
%--------------------------------------------------------------%

%--------------------------------------------------------------%
% define the thermal conductivity of the propagation medium
%--------------------------------------------------------------%
medium.thermal_conductivity = waterTc*ones(Nx, Ny);  % [m/s]

a = makeDisc(Nx, Ny, Nx/2, Ny/2, muscleRadius/dx);
medium.thermal_conductivity(a == 1) = muscleTc; 

b = 2*makeDisc(Nx, Ny, Nx/2, Ny/2, fatRadius/dx);
b(a == 1) = 1;
medium.thermal_conductivity(b == 2) = fatTc;

c = 3*makeDisc(Nx, Ny, Nx/2, Ny/2, skinRadius/dx);
c(a==1) = 1;
c(b==2) = 2;
medium.thermal_conductivity(c == 3) = skinTc;

f = 4*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, (boneRadius-boneThickness)/dx);
medium.thermal_conductivity(f == 4) = fatTc;

e = 5*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, boneRadius/dx);
e(f == 4) = 1;
medium.thermal_conductivity(e == 5) = boneTc;

d = makeDisc(Nx, Ny, Nx/2+tumorDistance/dx, Ny/2, tumorRadius/dx);
medium.thermal_conductivity(d == 1) = tumorTc;

medium.thermal_conductivity = repmat(medium.thermal_conductivity,[1,1,Nz]);
medium.thermal_conductivity(Nx/2+tumorDistance/dx-20:Nx/2+tumorDistance/dx+20,Ny/2-20:Ny/2+20,1:Nz/2-7) = muscleTc;
medium.thermal_conductivity(Nx/2+tumorDistance/dx-20:Nx/2+tumorDistance/dx+20,Ny/2-20:Ny/2+20,Nz/2+7:end) = muscleTc;
%define the perfusion rate of blood
medium.blood_perfusion_rate = zeros(Nx, Ny);
g = 6*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance1/dx, vesselRadius1/dx);  
medium.blood_perfusion_rate(g == 6) = bloodPr;
h = 7*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance2/dx, vesselRadius2/dx);  
medium.blood_perfusion_rate(h == 7) = bloodPr;
k = 8*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance3/dx, vesselRadius3/dx);  
medium.blood_perfusion_rate(k== 8) = bloodPr;
medium.blood_perfusion_rate = repmat(medium.blood_perfusion_rate,[1,1,Nz]);
%--------------------------------------------------------------%

%--------------------------------------------------------------%
% define the specific heat of the propagation medium
%--------------------------------------------------------------%
medium.specific_heat = waterSh*ones(Nx, Ny);  % [m/s]

a = makeDisc(Nx, Ny, Nx/2, Ny/2, muscleRadius/dx);
medium.specific_heat(a == 1) = muscleSh; 

b = 2*makeDisc(Nx, Ny, Nx/2, Ny/2, fatRadius/dx);
b(a == 1) = 1;
medium.specific_heat(b == 2) = fatSh;

c = 3*makeDisc(Nx, Ny, Nx/2, Ny/2, skinRadius/dx);
c(a==1) = 1;
c(b==2) = 2;
medium.specific_heat(c == 3) = skinSh;

f = 4*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, (boneRadius-boneThickness)/dx);
medium.specific_heat(f == 4) = fatSh;

e = 5*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, boneRadius/dx);
e(f == 4) = 1;
medium.specific_heat(e == 5) = boneSh;

d = makeDisc(Nx, Ny, Nx/2+tumorDistance/dx, Ny/2, tumorRadius/dx);
medium.specific_heat(d == 1) = tumorSh;

medium.specific_heat = repmat(medium.specific_heat,[1,1,Nz]);
medium.specific_heat(Nx/2+tumorDistance/dx-20:Nx/2+tumorDistance/dx+20,Ny/2-20:Ny/2+20,1:Nz/2-7) = muscleSh;
medium.specific_heat(Nx/2+tumorDistance/dx-20:Nx/2+tumorDistance/dx+20,Ny/2-20:Ny/2+20,Nz/2+7:end) = muscleSh;


medium.blood_specific_heat = zeros(Nx, Ny);
g = 6*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance1/dx, vesselRadius1/dx);  
medium.blood_specific_heat(g == 6) = bloodSh;
h = 7*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance2/dx, vesselRadius2/dx);  
medium.blood_specific_heat(h == 7) = bloodSh;
k = 8*makeDisc(Nx, Ny, Nx/2, Ny/2+vesselDistance3/dx, vesselRadius3/dx);  
medium.blood_specific_heat(k== 8) = bloodSh
medium.blood_specific_heat = repmat(medium.blood_specific_heat,[1,1,Nz]);
%--------------------------------------------------------------%

%%
% set the background temperature and heating terma
Q(:,:,1:Nz/2-tumorRadius/dx-1) = 0;
Q(:,:,Nz/2+tumorRadius/dx+1:end) = 0;
source.Q = 6.6*Q;   %17.5
T0 = 37*ones(Nx,Ny,Nz);
T0(medium.density==waterDensity) = 30;
source.T0 = T0;

% create binary sensor mask with a single sensor point in the centre of the
% grid
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1:Nx, 1:Ny, Nz/2) = 1;  

% create kWaveDiffusion object
input_args = {'PlotSim', false};
kdiff = kWaveDiffusion(kgrid, medium, source, sensor, input_args{:});

% set time step size
dt = 60;
on_time = 60;

% temperature controlling system 
min = 0;
while 1
    min = min + 1;
    kdiff.Q = 7.5*Q;
    kdiff.takeTimeStep(round(on_time / dt), dt);
    kdiff.T(medium.density == waterDensity) = 30;
    T_now = kdiff.T(round(Nx/2+tumorDistance/dx)-1,Ny/2,Nz/2);
    if(T_now >= 42.5)
        break
    end
end

% % =========================================================================
% % VISUALISATION
% % =========================================================================
% %
% % create time axis
t_axis = (0:kdiff.time_steps_taken - 1) / (60/dt);

% plot temperature time profile]
figure;
plot(t_axis, kdiff.sensor_data((Ny/2-1)*Ny+Nx/2+round(tumorDistance/dx)-1,:));
% xlabel('Time [min]');
% ylabel('Temperature [^\circC]');
% title('True-T');
filename = 'HeatingStage1.mat';save(filename);

clear;
% =========================================================================
% ACOUSTIC SIMULATION
% =========================================================================

% define the PML size
pml_size = 10;              % [grid points]

% define the grid parameters
Nx = 510 - 2 * pml_size;    % [grid points]
Ny = 510 - 2 * pml_size;    % [grid points]
dx = 220e-3/Nx;               % [m]
dy = 220e-3/Ny;               % [m]

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
muscleRadius = 6e-2;
fatRadius = muscleRadius + 1e-2;
skinRadius = fatRadius + 1.7e-3;
transRadius = 10e-2;
tumorRadius = 6.5e-3;
boneRadius = 1.25e-2;
boneThickness = 0.58e-2;

muscleSpeed = 1580;
waterSpeed = 1482;
fatSpeed = 1430;
marrowSpeed = 1371.9;
skinSpeed = 1595;
tumorSpeed = 1450;
boneSpeed = 2198;
bloodSpeed = 1570;
vesselWallSpeed = 1550;

boneDistance = 30e-3;
tumorDistance = 30e-3;

muscleAtten = 0.57;
fatAtten = 0.8; 
waterAtten = 2.17e-3;
skinAtten = 0.264;
boneAtten = 2.54;
tumorAtten = 0.4;
bloodAtten = 0.15;
vesselWallAtten = 0.6;
% tumorAtten = 0.75; %[dB/(MHz^y cm)]

muscleDensity = 1090;
waterDensity  = 1000;
fatDensity    = 911;
skinDensity   = 1109;
tumorDensity  = 1079;
boneDensity   = 1178;
bloodDensity = 1050;
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
% define the speed of the propagation medium
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
medium.sound_speed(f == 4) = fatSpeed+540;

e = 5*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, boneRadius/dx);
e(f == 4) = 1;
medium.sound_speed(e == 5) = boneSpeed;

d = makeDisc(Nx, Ny, Nx/2+tumorDistance/dx, Ny/2, tumorRadius/dx);
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
medium.alpha_coeff(f == 4) = fatAtten+1.54;
e = 5*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, boneRadius/dx);
e(f == 4) = 1;
medium.alpha_coeff(e == 5) = boneAtten;

d = makeDisc(Nx, Ny, Nx/2+tumorDistance/dx, Ny/2, tumorRadius/dx);
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
%--------------------------------------------------------------%
% define power 
%--------------------------------------------------------------%
medium.alpha_power = 1.5;
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
medium.density(f == 4) = fatDensity+120;
e = 5*makeDisc(Nx, Ny, Nx/2-boneDistance/dx, Ny/2, boneRadius/dx);
e(f == 4) = 1;
medium.density(e == 5) = boneDensity;

d = makeDisc(Nx, Ny, Nx/2+tumorDistance/dx, Ny/2, tumorRadius/dx);
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

% calculate the time step using an integer number of points per period
freq = 2e6;             % [Hz]
sampling_freq = 20e6;                        % period [s]
dt  = 1 / sampling_freq;                     % time step [s]

% calculate the number of time steps to reach steady state
sensor.record_start_index = 1;
Nt = 6000;

% create the time array
kgrid.setTime(Nt, dt);

%% define the source parameters(set transducer)
source.p_mask = zeros(Nx,Ny);
sensor.mask = zeros(Nx, Ny);
pieceAngle = 2*pi/256;
num1 = 2*pi/pieceAngle;
load('X_Pos.mat');
load('Y_Pos.mat');

for i = 1:num1
    sensor.mask(X_Pos(i), Y_Pos(i)) = 1;
end
% set the input arguements
input_args = {'PMLInside', false, 'PlotPML', false,'PlotSim', false, 'DisplayMask', ...
    'off'};

for i = 1:num1
    source.p0 = zeros(Nx, Ny);
    source.p0(X_Pos(i), Y_Pos(i)) = 10;
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(X_Pos(i), Y_Pos(i)) = 1;
%     sensor.mask = zeros(Nx, Ny);
%     sensor.mask(X_Pos(i), Y_Pos(i)) = 1;
    sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
    Datas{i} = sensor_data;
end
save('DataCollection_256_12101.mat','Datas');

%% Reconstruct B-mode image
% load('DataCollection_256_12101.mat');
B_mode = zeros(Nx, Ny);
m = 0;

for i = 1:num1
    tempB_mode1 = zeros(Nx, Ny);
    theta0Sin = (Ny/2 - Y_Pos(i))/(transRadius/dx);
    if(theta0Sin >= 0)
        theta0 = real(acos((Nx/2-X_Pos(i))/(transRadius/dx)));
    else
        theta0 = 2*pi - real(acos((Nx/2-X_Pos(i))/(transRadius/dx)));
    end
    theta1 = theta0 + pi/6;
    x1 = round(Nx/2 - (transRadius/dx)*cos(theta1));
    y1 = round(Ny/2 - (transRadius/dx)*sin(theta1));
    theta2 = theta0 - pi/6;
    x2 = round(Nx/2 - (transRadius/dx)*cos(theta2));
    y2 = round(Ny/2 - (transRadius/dx)*sin(theta2));
    xRangemin = min([x1 X_Pos(i) x2]);
    xRangemax = max([x1 X_Pos(i) x2]);
    yRangemin = min([y1 Y_Pos(i) y2]);
    yRangemax = max([y1 Y_Pos(i) y2]);
    iteration = 0;
    for n = 1:num1
        if(X_Pos(n)>=xRangemin && X_Pos(n)<=xRangemax && Y_Pos(n)>=yRangemin && Y_Pos(n)<=yRangemax)
            tempB_mode = zeros(Nx, Ny);
            iteration = iteration + 1;
            %106/524
            for j = 88:402
                for k = 88:402
                    %                         if(sqrt((j-Nx/2)^2 + (k-Ny/2)^2) < (transRadius/dx)-20)
                    %trans to point
                    routeLength1 = sqrt((X_Pos(i)-j)^2+(Y_Pos(i)-k)^2);
                    routeLength2 = sqrt((X_Pos(n)-j)^2+(Y_Pos(n)-k)^2);
                    routeLength = routeLength1 + routeLength2;
                    %                             num_dt = round(routeLength*dx/waterSpeed/dt) + m;
                    num_dt = routeLength*dx/1540/dt + m;
                    tempB_mode(j,k) = (Datas{i}(n,ceil(num_dt))-Datas{i}(n,floor(num_dt)))*num_dt+ceil(num_dt)*Datas{i}(n,floor(num_dt))-floor(num_dt)*Datas{i}(n,ceil(num_dt));
                    %                             tempB_mode(j,k) = Datas{i}(n,num_dt);
                    %                         tempB_mode(j,k) = interp1(round(num_dt-1):round(num_dt+1),Datas{i}(n,round(num_dt-1):round(num_dt+1)),num_dt);
                    %                         end
                end
            end
            tempB_mode1 = tempB_mode1 + tempB_mode;
        end
    end
    tempB_mode1 = tempB_mode1/iteration;
    B_mode = B_mode + tempB_mode1;
end
figure;imagesc(abs(B_mode));axis xy;
save('B_mode.mat','B_mode')
    
%% modify the image
v = B_mode;

for x = 1:490
    for y = 1:490
        v(x,y) = v(x,y) + 1.1;
        v(x,y) = log10(v(x,y));
        if((x-245)^2+(y-245)^2 >= 158^2)
            v(x,y) = 0;
        end
    end
end

figure;imagesc([0:22:220],[0:22:220],v);axis xy;
xlabel('x position (mm)');ylabel('y position (mm)');



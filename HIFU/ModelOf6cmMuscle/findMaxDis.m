clear;
load('Gridx_Pos.mat');
load('Gridy_Pos.mat');

Nx = 490;
Ny = 490;
dx = 220e-3/Nx;
tumorDistance = 25e-3;
% tumorDistance1 = 5e-3;
for i = 1:256

    Distance(i) = sqrt((x_Pos(i)-(Nx/2+tumorDistance/dx))^2 + (y_Pos(i)-(Ny/2))^2);

end
MaxDis = max(Distance);
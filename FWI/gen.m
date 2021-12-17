%% ========================================================================
% GENERAL INITIALIZATION
% =========================================================================
disk_name = 'E';
data_path = ':\FWI\data\';
load([disk_name,data_path,'medium200_blood.mat']);
medium200_blood(medium200_blood==1510)=1516.3;
load([disk_name,data_path,'sensor496_cord.mat']);

folderPath = [disk_name,':\FWI\gen_Data'];
cellSize = 1.3e-4;
numElements = 128;
sensor_freq = 20e6;
phySize = [0.02 0.02 ] * (3.125E-4/cellSize);
imgSize=[200;200];
%**************************************************************************
medium.sound_speed = medium200_blood;
medium.density = 1000 * ones(size(medium200_blood));

%% ========================================================================
sensorX = sourceXrad_496(1:3.87:end);
sensorY = sourceYrad_496(1:3.87:end);
sensormask =zeros(size(medium200_blood));
for i = 1:size(sensorX,1)
   sensormask(sensorX(i),sensorY(i)) =1;
end
[sensorX,sensorY]=find(sensormask~=0);
% SIMULATION OF RECORDED WAVE FIELDS
[waveDataExp, gradData, tArr, elemArr] = waveSim_ring_CPU_FWI(numElements, medium, cellSize, imgSize, 0, 0, 0,...
    sensorX,sensorY,sensor_freq);

% save the wavefield input
%********************************************************************
save(sprintf('%s/%s', folderPath, 'forwarddata_01'), 'waveDataExp', 'tArr', 'phySize', 'imgSize', 'medium',...
    'numElements','cellSize', 'imgSize');
%========================================================

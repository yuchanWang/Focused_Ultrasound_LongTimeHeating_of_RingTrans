
%% ========================================================================
% GENERAL INITIALIZATION
% =========================================================================
%**********************************************************
disk_name = 'E';
data_path = ':\FWI\data\';
gen_path = ':\FWI\gen_Data\';
load([disk_name,data_path,'medium200_blood.mat']);
load([disk_name,data_path,'sensor496_cord.mat']);
load([disk_name,data_path,'roi_blood.mat']);
%**********************************************************
load([disk_name,gen_path,'forwarddata_01.mat']);
sensorX = sourceXrad_496(1:3.87:end);
sensorY = sourceYrad_496(1:3.87:end);
sensormask =zeros(200,200);
for i = 1:size(sensorX,1)
   sensormask(sensorX(i),sensorY(i)) =1;
end
[sensorX,sensorY]=find(sensormask~=0);

% make directory to store reconstructed images
savedPath = [disk_name,':\FWI\Recons\'];
%**********************************************************
subfolder = sprintf('recdata');
sensor_freq = 20e6;

mkdir([savedPath, subfolder]);
%********************************************************
roi = makeDisc(200, 200, 100, 100, 58);
roi(roi_blood==5)=0;
roi(roi_blood==6)=0;
roi(roi_blood==7)=0;

roiMap = roi_blood;
updatedIdx = logical(roi(:));

% create the density map (kg/m^3)
medium.density = 1000 * ones(size(roiMap)); % 

% other optimization parameters
TOL = 0.01;
%% ========================================================================
Cd = 0.001^2;
Cdinv = 1/Cd;
% covariance matrix for model matrix and correlation adjustments
numRoi = max(roiMap(:)); 
corrCoeff = zeros(1, numRoi);
sigmaRoi = zeros(1, numRoi); 
corrCoeff(1:numRoi) = 0.0005;
sigmaRoi(1) = 10;
sigmaRoi(2) = 10;
sigmaRoi(3) = 10;
sigmaRoi(4) = 10;
sigmaRoi(5) = 10;
sigmaRoi(6) = 10;
sigmaRoi(7) = 10;
sigmaRoi(8) = 10;
%% ========================================================================
% INITIALIZATION FOR 1st ITERATION
% =========================================================================
%*******************************************
initSOSMap = 1450 * ones(size(medium200_blood));
initSOSMap(medium200_blood==1482)=1482; 
initSOSMap(roi_blood==6)=1371.9; 
initSOSMap(medium200_blood==1550)=1550; 
initSOSMap(medium200_blood==2198)=2198; 
initSOSMap = initSOSMap./1000;
priorInfo.map = initSOSMap * 1000;
medium.sound_speed = priorInfo.map;
SpeedImage = medium.sound_speed;

%% ========================================================================
% ITERATIVE RECONSTRUCTION
% =========================================================================
uk = SpeedImage(:)/1000;
uk = uk(updatedIdx);
mk = SpeedImage(:)/1000;
mk = mk(updatedIdx);
mPrior = priorInfo.map(updatedIdx)/1000;

i = 0;
numItr = 301;
errorArr = [];
%***********************************
stepsize = 1;

while (i < numItr)
    if(i > 0 && resumeSim == 0)
        epsilon = max(abs(mk))/max(abs(gk))/500*stepsize;
        if(epsilon > (1*stepsize)) 
            epsilon = 1*stepsize;
        end
        
        fprintf('\n2) Performing line search...\n');

        mkInc = mk + epsilon*gk';
        medium.sound_speed(updatedIdx) = mkInc*1000;
        medium.sound_speed(medium.sound_speed <= 5) = 5;
        [waveDataIncSim, gradData, tArr, elemArr] = waveSim_ring_CPU_FWI(numElements, medium, cellSize, imgSize, 0, 0, 0,...
            sensorX,sensorY,sensor_freq);
        Jkdk = (waveDataIncSim(:)-waveDataSim(:))/epsilon;
        
        a = gk * invcovmv(roiMap(updatedIdx), J', 1:numRoi, corrCoeff, sigmaRoi);
        c = gk * invcovmv(roiMap(updatedIdx), gk', 1:numRoi, corrCoeff, sigmaRoi);
        b = Jkdk' * Cdinv * Jkdk; 
        betak = -a / (b+c);
        
        mkPrior = mk;
        mk = mk + betak*gk';
        medium.sound_speed(updatedIdx) = mk*1000;
        SpeedImage = medium.sound_speed;
        
        fprintf('\n3) Saving results...\n');
        
        save(sprintf('%s%s/iteration_%i', savedPath, subfolder, i), 'SpeedImage', 'J', 'Jd', 'gk', 'betak');
        
        figure(1);
        plot(errorArr);
        title('Gradient magnitude of step  (L2-COV)');
        drawnow;
    end
    
    fprintf('\n=== Iteration %i ===\n', i + 1);
    
    % forward and adjoint solvers on speed of sound image, given a pulse as a source point
    fprintf('\n1) Solving forward model and calculating gradient...\n');
    [waveDataSim, Jd, tArr, elemArr] = waveSim_ring_CPU_FWI(numElements, medium, cellSize, imgSize, 1, 0, waveDataExp,...
        sensorX,sensorY,sensor_freq);
    Jd = Jd(updatedIdx);
    Jd = Jd*1000;
    if(i == 0)
        J =  Cdinv * covmv(roiMap(updatedIdx), Jd, 1:numRoi, corrCoeff, sigmaRoi) + (mk - mPrior)';
        gk = -J;
    else
        JPrior = J;
        
        J =  Cdinv * covmv(roiMap(updatedIdx), Jd, 1:numRoi, corrCoeff, sigmaRoi) + (mk - mPrior)';
        
        gamma1 = invcovmv(roiMap(updatedIdx), J - JPrior, 1:numRoi, corrCoeff, sigmaRoi) ...
                        * J';
        gamma2 = invcovmv(roiMap(updatedIdx), JPrior, 1:numRoi, corrCoeff, sigmaRoi) ...
                        * JPrior';

        gamma = max(0, gamma1/gamma2);
        
        gk = -J + gamma*gk;
    end

    errVal = sqrt(sum(J.^2));
    i = i + 1;

    if(i == 1)
        e0 = errVal;
    end
    
    errorArr = [errorArr errVal]; 
    
    resumeSim = 0;
end
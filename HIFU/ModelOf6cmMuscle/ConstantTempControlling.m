clear;
load('HeatingStage1.mat');
deltak1 = 0.01;
deltak2 = 0.01;
for i = (min+1):(300+min)      %9:188/13:192/
    T_last2 = kdiff.sensor_data((Ny/2-1)*Ny+Nx/2+round(tumorDistance/dx)-1,i-2);
    T_last1 = kdiff.sensor_data((Ny/2-1)*Ny+Nx/2+round(tumorDistance/dx)-1,i-1);
    
    Q_last = kdiff.Q(Nx/2+round(tumorDistance/dx)-1,Ny/2,Nz/2);
    
    Q_source = kdiff.Q;

    deltaT = abs(T_last1-T_last2);
    
    if(deltaT < 0.001)
        k = -7.5e6;
    end
    if((deltaT < 0.01) && (deltaT >= 0.001))
        k = -3e6;  %-2.5e6

    end
    if((deltaT < 0.1) && (deltaT >= 0.01))
        k = -7.4e5;  %-8.5e5´ó
    end
    if(deltaT >= 0.1)
        k = -1e5; %-6e5´ó
    end
    
    Q_new = k*(T_last1 - T_last2) + Q_last;
    if(T_last1 > 43.05)  %42.57&42.45
        k = k - deltak1*k;   %0.03
        deltak1 = deltak1 + 0.01;
        Q_new = k*(T_last1 - T_last2) + Q_last - 0.0015*Q_last;
    end
    
    if(T_last1 < 42.95 && i > 30)  %42.57&42.45
        k = k + deltak2*k;   %0.03
        deltak2 = deltak2 + 0.01;
        Q_new = k*(T_last1 - T_last2) + Q_last + 0.001*Q_last;
    end
    
    kdiff.Q = (Q_new/Q_last).*Q_source;
    kdiff.takeTimeStep(round(on_time / dt),dt);
    kdiff.T(medium.density == waterDensity) = 30;
end

% =========================================================================
% VISUALISATION
% =========================================================================
%%
% create time axis
t_axis = (0:kdiff.time_steps_taken - 1) / (60/dt);

% plot temperature time profile]
figure;
plot(t_axis, kdiff.sensor_data((Ny/2-1)*Ny+Nx/2+round(tumorDistance/dx)-1,:));hold on
plot([0,319],[43.20,43.20],'--r');hold on
plot([0,319],[42.80,42.80],'--r');hold on

% plot(t_axis, kdiff.sensor_data);
xlabel('Time [min]');
ylabel('Temperature [^\circC]');
% title('True-T');

% figure;
% plot(t_axis, Fault_temp);
% xlabel('Time [min]');
% ylabel('Temperature [^\circC]');
% title('False-T');


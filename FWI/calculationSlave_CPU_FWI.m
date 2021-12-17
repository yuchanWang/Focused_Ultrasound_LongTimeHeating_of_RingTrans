function [sensor_data, grad_data, norm_data] = calculationSlave_CPU_FWI(kgrid, medium, PMLSize, source_sim, sensor, ...
    input_args, rcvGridIdx, adjoint, waveDataExp, dt, elNum)
  
    fprintf('On element %i... ', elNum);
    [sensor_data] = kspaceFirstOrder2D(kgrid, medium, source_sim, sensor, input_args{:});
	sensor_data = sensor_data';
    
    waveDiff = sensor_data(:, rcvGridIdx) - waveDataExp;
    waveDiff = waveDiff(end:-1:1, :);
    
    if(adjoint) % also compute the adjoint data
        sensor.mask = ones(kgrid.Nx, kgrid.Ny);
        
        source_adj.p_mask = zeros(kgrid.Nx, kgrid.Ny);

        source_adj.p_mode = 'dirichlet'; % enforce dirichlet boundary condition
        sourceGridIdx = rcvGridIdx; % highlight the region that all sources are in
        source_adj.p_mask(sourceGridIdx) = 1;
        source_adj.p = zeros(length(rcvGridIdx), size(waveDiff(:, :), 1));

        % for each element, compute adjoint field
        for i = 1:length(rcvGridIdx)
            source_adj.p(i, :) = waveDiff(:, i)/4/pi; % set the residual wavefields to all receivers
        end
        
        input_args = {'DisplayMask', source_adj.p_mask, 'PMLInside', false, 'PlotPML', false, ...
        'PMLSize', PMLSize, 'PlotSim', false, 'RecordMovie', false, 'Smooth', [true, true, true], ...
        'DataCast', 'single'};

        [out_data] = kspaceFirstOrder2D(kgrid, medium, source_adj, sensor, input_args{:});
        out_data = out_data';

        Jd = zeros(1, size(medium.sound_speed, 1) * size(medium.sound_speed, 2));
        L = size(sensor_data, 1);
        for l = 2:L-1
            Jd = Jd + 2./(medium.sound_speed(:).^3)' .* out_data(L - l + 1, :) .* ...
                (sensor_data(l - 1, :) - 2 * sensor_data(l, :) + sensor_data(l + 1, :))/dt;
        end
        
        grad_data = Jd;
        
        norm_data = zeros(1, size(medium.sound_speed, 1) * size(medium.sound_speed, 2));
        
        for l = 1:L
            norm_data = norm_data + (sensor_data(l, :).^2)*dt;
        end
    else
        grad_data = zeros(1, size(medium.sound_speed, 1) * size(medium.sound_speed, 2));
        norm_data = zeros(1, size(medium.sound_speed, 1) * size(medium.sound_speed, 2));
    end
    
    sensor_data = sensor_data(:, rcvGridIdx);
end
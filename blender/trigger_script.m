clear all

load('sar_consts.mat');

model = 'sar_model';
load_system(model);
simMode = get_param(model, 'SimulationMode');
flightTime = str2double(get_param(model, 'StopTime'));

reference_range = 800;

sim_out = sim(model, 'SimulationMode', simMode);

figure(1)
imagesc(real(reshape(sim_out.spotlight_data, [], 4001))); title('SAR Raw Spotlight Mode Data')
xlabel('Cross-range')
ylabel('Down-range')

figure(2)
imagesc(real(reshape(sim_out.stripmap_data, [], 4001))); title('SAR Raw Stripmap Mode Data')
xlabel('Cross-range')
ylabel('Down-range')

matched_filter_coeffs = sim_out.lin_fm_coeffs(1:360);

pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', 'PropagationSpeed', c, 'SampleRate', fs);
[spotlight_data, spotlight_range_grid] = pulseCompression(sim_out.spotlight_data, matched_filter_coeffs);
[stripmap_data, stripmap_range_grid] = pulseCompression(sim_out.stripmap_data, matched_filter_coeffs);

spotlight_reshaped_data = reshape(spotlight_data, 2002, []);
reshaped_spotlight_range_grid = reshape(spotlight_range_grid, 2002, []);
stripmap_reshaped_data = reshape(spotlight_data, 2002, []);
reshaped_stripmap_range_grid = reshape(spotlight_range_grid, 2002, []);

figure(3)
imagesc(real(spotlight_reshaped_data)); title('SAR Range Compressed Data - Spotlight')
xlabel('Cross-range')
ylabel('Down-range')

figure(4)
imagesc(real(stripmap_reshaped_data)); title('SAR Range Compressed Data - Stripmap')
xlabel('Cross-range')
ylabel('Down-range')

omega_k_image_spotlight = omegak(spotlight_reshaped_data, fs, maxRange, fc, flightTime, speed, reference_range, prf);
omega_k_image_stripmap = omegak(stripmap_reshaped_data, fs, maxRange, fc, flightTime, speed, reference_range, prf);

figure(5)
imagesc((abs(omega_k_image_spotlight.')));
title('Spotlight SAR Data focused using Omega-K algorithm')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

figure(6)
imagesc((abs(omega_k_image_stripmap.')));
title('Stripmap SAR Data focused using Omega-K algorithm')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

backproject_image_spotlight = helperBackProjection(spotlight_reshaped_data, reshaped_spotlight_range_grid, fastTime, fc, fs, speed, crossRangeResolution, c, prf);
backproject_image_stripmap = helperBackProjection(stripmap_reshaped_data, reshaped_stripmap_range_grid, fastTime, fc, fs, speed, crossRangeResolution, c, prf);

figure(7)
imagesc((abs(backproject_image_spotlight)))
title('Spotlight SAR Data focused using back projection algorithm')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

figure(8)
imagesc((abs(backproject_image_stripmap)))
title('Spotlight SAR Data focused using back projection algorithm')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

function data = helperBackProjection(sigdata,rnggrid,fastTime,fc,fs,speed,crossRangeResolution,c, prf)
    data = zeros(size(sigdata));
    azimuthDist = -200:speed/prf:200;%azimuth distance
    rangelims = [700 1400];
    crossrangelims = [-10 10];
    rangeIdx =  [find(rnggrid>rangelims(1), 1) find(rnggrid<rangelims(2),1,'last')];
    crossrangeIdxStart = find(azimuthDist>crossrangelims(1),1);
    crossrangeIdxStop = find(azimuthDist<crossrangelims(2),1,'last');
    for i= rangeIdx(1):rangeIdx(2)  % Iterate over the range indices
        
        % Using desired cross-range resolution, compute the synthetic aperture
        % length
        lsynth= (c/fc)* (c*fastTime(i)/2)/(2*crossRangeResolution);
        lsar = round(lsynth*length(azimuthDist)/azimuthDist(end)) ;
        % Ensure lsar is an odd number
        lsar = lsar + mod(lsar,2);    
        
        % Construct hanning window for cross-range processing, to suppress the
        % azimuthal side lobes
        hn= hanning(lsar).';
        % Iterate over the cross-range indices
        for j= crossrangeIdxStart:crossrangeIdxStop 
            % azimuth distance in x direction over cross-range indices
            posx= azimuthDist(j);
            % range in y-direction over range indices
            posy= c*fastTime(i)/2;
            % initializing count to zero
            count= 0;
            % Iterate over the synthetic aperture
            for k= j-lsar/2 +1:j+ lsar/2 
                % Time delay for each of range and cross-range indices
                k_round = round(k/100);
                td= sqrt((azimuthDist(k)- posx)^2 + posy^2)*2/c;
                cell= round(td*fs) +1 ;
                signal = sigdata(cell, k_round);
                count= count + hn(k -(j-lsar/2))*signal *exp(1j*2*pi*fc*(td));
            end
            % Processed data at each of range and cross-range indices
            data(i,j)= count;
        end
        
    end
end
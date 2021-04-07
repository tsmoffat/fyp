clear all

load('sar_consts.mat');

model = 'sar_model';
load_system(model);
simMode = get_param(model, 'SimulationMode');
flightTime = str2double(get_param(model, 'StopTime'));

reference_range = 800;

sim_out = sim(model, 'SimulationMode', simMode);

figure(3)
imagesc(real(reshape(sim_out.range_data, [], 4001))); title('SAR Raw Data')
xlabel('Cross-range')
ylabel('Down-range')

matched_filter_coeffs = sim_out.lin_fm_coeffs(1:360);

pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', 'PropagationSpeed', c, 'SampleRate', fs);
[cdata, range_grid] = pulseCompression(sim_out.range_data, matched_filter_coeffs);

reshaped_data = reshape(cdata, 2002, []);
reshaped_range_grid = reshape(range_grid, 2002, []);

figure(4)
imagesc(real(reshaped_data)); title('SAR Range Compressed Data')
xlabel('Cross-range')
ylabel('Down-range')

omega_k_image = omegak(reshaped_data, fs, maxRange, fc, flightTime, speed, reference_range, prf);

figure(1)
imagesc((abs(omega_k_image.')));
title('SAR Data focused using Omega-K algorithm')
xlabel('Cross-Range Samples')
ylabel('Range Samples')

back_project_image = helperBackProjection(reshaped_data, reshaped_range_grid, fastTime, fc, fs, speed, crossRangeResolution, c, prf);

figure(2)
imagesc((abs(back_project_image)))
title('SAR Data focused using back projection algorithm')
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
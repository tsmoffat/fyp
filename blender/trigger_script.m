clear all

load('sar_consts.mat');

model = 'sar_model';
load_system(model);
simMode = get_param(model, 'SimulationMode');
flightTime = str2double(get_param(model, 'StopTime'));

reference_range = 800;
cross_range_resolution = 1;

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

backproject_image_spotlight = backprojection(spotlight_reshaped_data, reshaped_spotlight_range_grid, fastTime, prf, speed, cross_range_resolution, fc, fs);
backproject_image_stripmap = backprojection(spotlight_reshaped_data, reshaped_spotlight_range_grid, fastTime, prf, speed, cross_range_resolution, fc, fs);

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
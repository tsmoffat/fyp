function omega_k_output = omegak(raw_data, fs, maxRange, fc, flightTime, speed, reference_range, prf)
    c = physconst('LightSpeed');
    truncated_range_samples = ceil((2*maxRange/c)*fs);
    pulses = length(0:1/fs:(truncated_range_samples-1)/fs);

    frequency_range = linspace(fc-fs/2, fc + fs/2, pulses);
    range_delta = 2*(2*pi*frequency_range)/c;
    cross_range_wavenumber = 2*pi*linspace(-prf/2, prf/2, (flightTime/(1/prf) + 1))./speed;

    cross_range_number_matrix = cross_range_wavenumber.';

    temp_wavenumber = range_delta.^2-cross_range_number_matrix.^2;
    temp_wavenumber = sqrt(temp_wavenumber.*(temp_wavenumber>0));
    final_wavenumbers_for_focusing = exp(1i*temp_wavenumber.*reference_range);
    fft_data =fftshift(fft(fftshift(fft(raw_data,[],1),1),[],2),2);
    %fft_data = fftshift(fft2(raw_data));

    bulk_compression = (fft_data.').*final_wavenumbers_for_focusing;

    stolt_interpolation = bulk_compression;
    for i = 1:size((bulk_compression),1)
        stolt_interpolation(i,:) = interp1(temp_wavenumber(i,:), bulk_compression(i,:), range_delta(1,:));
    end

    stolt_interpolation(isnan(stolt_interpolation)) = 1e-30;

    stolt_interpolation = stolt_interpolation.*exp(-1i*range_delta.*reference_range);

    omega_k_output = ifft2(stolt_interpolation);


end
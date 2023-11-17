function TT_synth = create_synthetic_dataset(datetime_vec, wave_starts, wave_ends, wave_periods, wave_amplitudes, noise_floor)
%CREATE_SYNTHETIC_DATASET Summary of this function goes here
%   Detailed explanation goes here
rng(1)

baseline_pres = 1000 + 2*randn;
baseline_pres_vec = repmat(baseline_pres, size(datetime_vec));

noise_vec = noise_floor * randn(size(datetime_vec));

full_pres_vec = baseline_pres_vec + noise_vec;

if isduration(wave_periods)
    wave_periods = seconds(wave_periods);
end

for cur_wave = 1:length(wave_starts)
    wave_start_ind = find(datetime_vec == wave_starts(cur_wave));
    wave_end_ind = find(datetime_vec == wave_ends(cur_wave));

    wave_x = (0:(wave_end_ind-wave_start_ind)) .* (2*pi/wave_periods(cur_wave));

    % ramp up the wave amplitude until 10% of the event duration, then hold
    % constant for 80% of the event duration, then ramp down for the last
    % 10% of the event duration
    rampup = linspace(0, 1, ceil(length(wave_x)/10));
    rampdown = linspace(1, 0, ceil(length(wave_x)/10));
    plateau = ones(1, length(wave_x) - length(rampup) - length(rampdown));
    amplitude_curve = [rampup, plateau, rampdown];

    wave_trace =  wave_amplitudes(cur_wave) * amplitude_curve .* sin(wave_x);

    full_pres_vec(wave_start_ind:wave_end_ind) = full_pres_vec(wave_start_ind:wave_end_ind) + wave_trace;

end
synth_Temp = repmat(20, size(full_pres_vec));

TT_synth = timetable(datetime_vec(:), synth_Temp(:), full_pres_vec(:));

end


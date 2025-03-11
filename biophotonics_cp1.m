%%% Biophotonics Checkpoint 1 Data Analysis:
%%% Take note: there are body tremors at: 3,5, 9, 13Hz, as well as 
%%% at very low frequencies (low pass filter)->band pass dc noise

% Part 1: Load the Signal Data from the Oscilloscope
file = "randomname.dat";
data = load(file);
signal = data(:, 2);                                                        % Second .dat data is the signal amplitude
t = data(:, 1);                                                             % First .dat data is the time passed

% Part 2: Compute the Sampling Frequency and # of Samples
fs = 1 / mean(diff(t));
N = length(signal);

% Part 3: Perform FFT on the Signal
Y = fft(signal);
freq = (0 : N-1) * (fs / N);
Y_magnitude = abs(Y(1 : N/2));
freq_pos = f(1 : N/2);
r
Y_spec = Y_magnitude ./ freq_pos;

% Part 4: Remove DC Noise and Apply Bandpass Filter 
bp_low_freq = 0.5;
bp_high_freq = 5;

bp_signal = bandpass(signal, [bp_low_freq, bp_high_freq], fs);

% Part 5: Apply Narrowband Filter at 1Hz
focus_freq = 1;
freq_dev = 0.5;

nb_signal = bandpass(bp_signal, [focus_freq - freq_dev, focus_freq + freq_dev], fs);

% Part 5.b: Filter out any body tremor frequencies
bt_1 = 3;
nb_signal1 = bandpass(nb_signal, [bt_1 - freq_dev, bt_1 + freq_dev], fs);

bt_2 = 5;
nb_signal2 = bandpass(nb_signal1, [bt_2 - freq_dev, bt_2 + freq_dev], fs);

bt_3 = 9;
nb_signal3 = bandpass(nb_signal2, [bt_3 - freq_dev, bt_3 + freq_dev], fs);

bt_4 = 13;
nb_signal4 = bandpass(nb_signal3, [bt_4 - freq_dev, bt_4 + freq_dev], fs);


% Part 6: Plot the Original Oscilloscope vs. Normalized Spectrum vs.
% Narrowband Filtering
figure;
subplot(3, 1, 1);
plot(t, signal);
title("Original Oscilloscope Signal: ");
xlabel("Time (s)");
ylabel("Amplitude");

subplot(3, 1, 2);
plot(f_pos, Y_norm);
title("Normalized Spectrum: ");
xlabel("Frequency (Hz)");
ylabel("Magnitude");

subplot(3, 1, 2);
plot(t, nb_signal4);
title("Narrowband Filtered Signal at 1Hz: ");
xlabel("Time (s)");
ylabel("Amplitude");































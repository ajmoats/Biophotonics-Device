%% MATLAB Implementation: NIRS Data Processing from Oscilloscope
% This script processes two oscilloscope data files (630nm and 930nm)
% to compute the Metabolic Index (MI) for non-invasive blood glucose estimation.

clc; clear; close all;

%% Load Oscilloscope Data
file_630nm = 'ocp2.dat'; % Update with actual filename
file_930nm = 'ocp3.dat'; % Update with actual filename

% Load data assuming two columns: [Time, Signal]
data_630 = readmatrix(file_630nm);
data_930 = readmatrix(file_930nm);

% Extract time and signal values
time_630 = data_630(:,1);
signal_630 = data_630(:,2);

time_930 = data_930(:,1);
signal_930 = data_930(:,2);

% Ensure both signals have the same time vector
fs = 1 / (time_630(2) - time_630(1)); % Sampling frequency
if length(time_630) ~= length(time_930)
    error('Data lengths for 630nm and 930nm do not match. Ensure synchronization.');
end

%% Band-Pass Filtering to Extract AC Components
low_cutoff = 0.8;  % Hz (Heart rate lower bound)
high_cutoff = 10;  % Hz (Remove high-frequency noise)
[b, a] = butter(2, [low_cutoff high_cutoff] / (fs / 2), 'bandpass');

signal_630_AC = filtfilt(b, a, signal_630);
signal_930_AC = filtfilt(b, a, signal_930);

%% FFT and Phase Delay Calculation
y1 = fft(signal_630_AC);
y2 = fft(signal_930_AC);

N = length(signal_630);
freqs = (0:N-1) * (fs / N);

% Find dominant frequency (heart rate peak)
[~, peak_idx] = max(abs(y1));
peak_freq = freqs(peak_idx);

% Compute Phase Difference
phase_630 = angle(y1(peak_idx));
phase_930 = angle(y2(peak_idx));
delta_theta = phase_930 - phase_630;

%% Compute SaO2 and Metabolic Index (MI)
A_630 = max(signal_630_AC) - min(signal_630_AC);
A_930 = max(signal_930_AC) - min(signal_930_AC);
SaO2 = A_630 / (A_630 + A_930);

MI = SaO2 * (1 - SaO2) * abs(delta_theta);

%% Part 5: Apply Narrowband Filter at 1Hz
focus_freq = 1;
freq_dev = 0.5;

nb_signal = bandpass(signal_630_AC, [focus_freq - freq_dev, focus_freq + freq_dev], fs);

% Part 5.b: Filter out any body tremor frequencies
bt_1 = 3;
nb_signal1 = bandpass(nb_signal, [bt_1 - freq_dev, bt_1 + freq_dev], fs);

bt_2 = 5;
nb_signal2 = bandpass(nb_signal1, [bt_2 - freq_dev, bt_2 + freq_dev], fs);

bt_3 = 9;
nb_signal3 = bandpass(nb_signal2, [bt_3 - freq_dev, bt_3 + freq_dev], fs);

bt_4 = 13;
nb_signal4 = bandpass(nb_signal3, [bt_4 - freq_dev, bt_4 + freq_dev], fs);

%% Plot Results
figure;
subplot(2,1,1);
plot(time_630, signal_630_AC, 'r', time_930, signal_930_AC, 'b');
legend('630nm (Oxyhemoglobin)', '930nm (Deoxyhemoglobin)');
title('Filtered NIRS Signals'); xlabel('Time (s)'); ylabel('Amplitude');

subplot(2,1,2);
plot(peak_freq, abs(delta_theta), 'ko');
title('Phase Delay Analysis'); xlabel('Frequency (Hz)'); ylabel('Phase Delay (radians)');

fprintf('Computed Metabolic Index (MI): %.4f\n', MI);

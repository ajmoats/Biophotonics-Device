%%% Biophotonics Checkpoint X Data Analysis:
%%% Objective: Clean and Analyze Downward Crossings of Biophotonic Signal
%%% Note: Body tremors at 3, 5, 9, 13 Hz, low frequency drift -> low-pass recommended

% Red Signal
file = "C:\Users\keyaa\Downloads\adc_red_LJA.csv";
data = readmatrix(file);

x = data(:, 1); % better: time vector in seconds
y = data(:, 2); % better: raw amplitude signal (red channel)

% Known sampling frequency
Fs = 50; % Hz (based on your acquisition)
fs = 1 / mean(diff(x)); 
disp(['Using Sampling Frequency: ', num2str(Fs), ' Hz']); % better: confirm actual spacing vs known

% Raw Plot for Reference
figure;
plot(x, y, 'k', 'DisplayName', 'Raw Signal');
title("Raw Signal Plot, Red"); % better: visualize unprocessed red signal

% Design a Bandpass Filter
low_freq = 0.2;
high_freq = 3;
order = 1;
Wn = [low_freq high_freq] / (fs / 2);
[b, a] = butter(order, Wn, 'Bandpass');
y_filtered = filtfilt(b, a, y); % better: apply zero-phase bandpass filter

% Raw Bandpass for Reference
figure;
plot(x, y_filtered, 'r', 'DisplayName', 'Bandpass Signal');
title("Bandpass Signal Plot, Red"); % better: visualize cleaned red signal

% Part 3: Do Peak-Trough Detection for Reference
% Peak and Valley Detection
min_peak_dist = round(0.1 * Fs); % better: enforce minimum 100ms spacing between beats
prominence = 0.005; % better: ignore small fluctuations

[peak_vals, peak_locs] = findpeaks(y_filtered, 'MinPeakDistance', min_peak_dist, 'MinPeakProminence', prominence);
[valley_vals, valley_locs] = findpeaks(-y_filtered, 'MinPeakDistance', min_peak_dist, 'MinPeakProminence', prominence);
valley_vals = -valley_vals; % better: correct valley inversion

% Plot Peaks and Valleys with Intervals
figure;
plot(x, y_filtered, 'b');
hold on;
plot(x(peak_locs), peak_vals, 'ko', 'DisplayName', 'Peaks');
plot(x(valley_locs), valley_vals, 'go', 'DisplayName', 'Valleys');
for i = 2:length(peak_locs)
    t0 = x(peak_locs(i-1));
    t1 = x(peak_locs(i));
    interval = t1 - t0;
    midpoint = (t0 + t1) / 2;
    text(midpoint, peak_vals(i)+0.02, sprintf('%.2fs', interval), 'HorizontalAlignment', 'center', 'FontSize', 8);
end
legend;
title(sprintf('Detected Peaks and Valleys (Prominence = %.4f), Red', prominence));
xlabel('Time (s)');
ylabel('Filtered Amplitude (V)');
grid on;
title("Interval Signal Plot, Red");

% Infrared (IR) Signal
ir_file = "C:\Users\keyaa\Downloads\adc_ir_LJA.csv";
ir_data = readmatrix(ir_file);

ir_x = ir_data(:, 1); % better: IR time vector (seconds)
ir_y = ir_data(:, 2); % better: IR raw signal amplitude

% Known sampling frequency
Fs = 50; % Hz (based on your acquisition)
disp(['Using Sampling Frequency: ', num2str(Fs), ' Hz']);

% Raw Plot for Reference
figure;
plot(ir_x, ir_y, 'k', 'DisplayName', 'Raw Signal');
title("Raw Signal Plot, IR");

% Design a Bandpass Filter
ir_low_freq = 0.2;
ir_high_freq = 3;
ir_order = 1;
ir_Wn = [ir_low_freq ir_high_freq] / (fs / 2);
[ir_b, ir_a] = butter(ir_order, ir_Wn, 'Bandpass');
ir_y_filtered = filtfilt(ir_b, ir_a, ir_y); % better: apply bandpass to IR

% Raw Bandpass for Reference
figure;
plot(ir_x, ir_y_filtered, 'r', 'DisplayName', 'Bandpass Signal');
title("Bandpass Signal Plot, IR");

% Part 3: Do Peak-Trough Detection for Reference
% Peak and Valley Detection
ir_min_peak_dist = round(0.1 * Fs); % better: enforce minimum beat gap
ir_prominence = 0.005;

[ir_peak_vals, ir_peak_locs] = findpeaks(ir_y_filtered, 'MinPeakDistance', ir_min_peak_dist, 'MinPeakProminence', ir_prominence);
[ir_valley_vals, ir_valley_locs] = findpeaks(-ir_y_filtered, 'MinPeakDistance', ir_min_peak_dist, 'MinPeakProminence', ir_prominence);
ir_valley_vals = -ir_valley_vals;

% Plot Peaks and Valleys with Intervals
figure;
plot(ir_x, ir_y_filtered, 'b');
hold on;
plot(ir_x(ir_peak_locs), ir_peak_vals, 'ko', 'DisplayName', 'Peaks');
plot(ir_x(ir_valley_locs), ir_valley_vals, 'go', 'DisplayName', 'Valleys');
for ir_i = 2:length(ir_peak_locs)
    ir_t0 = ir_x(ir_peak_locs(ir_i-1));
    ir_t1 = ir_x(ir_peak_locs(ir_i));
    ir_interval = ir_t1 - ir_t0;
    ir_midpoint = (ir_t0 + ir_t1) / 2;
    text(ir_midpoint, ir_peak_vals(ir_i)+0.02, sprintf('%.2fs', ir_interval), 'HorizontalAlignment', 'center', 'FontSize', 8);
end
legend;
title(sprintf('Detected Peaks and Valleys (Prominence = %.4f)', ir_prominence));
xlabel('Time (s)');
ylabel('Filtered Amplitude (V)');
grid on;
title("Interval Signal Plot, IR");

% Attempt to plot both signals together
figure;
plot(x, y_filtered, 'k', 'DisplayName', 'Red Bandpass Signal');
hold on;
plot(ir_x, ir_y_filtered, 'r', 'DisplayName', 'IR Bandpass Signal');
title("Combined Red & IR Bandpassed Signal Plot")
hold off;

% Find the peaks in both samples and try to align them
% Step 1: Find Peaks
[~, locs_red] = findpeaks(y_filtered, 'MinPeakProminence', 0.01);
[~, locs_ir] = findpeaks(ir_y_filtered, 'MinPeakProminence', 0.01);

% Step 2: Convert Indices to Time
t_red_peaks = x(locs_red);
t_ir_peaks = ir_x(locs_ir);

% Step 3: Match and Calculate Delays
N = min(length(t_red_peaks), length(t_ir_peaks)); % better: truncate to matching peak count
delays = t_ir_peaks(1:N) - t_red_peaks(1:N); % better: IR minus Red delay vector
avg_delay = mean(delays); % better: average delay in seconds
disp(['Average Time Delay (IR - Red): ', num2str(avg_delay), ' seconds']);

% Step 4: Replot Signals with Peak Markers
figure;
plot(ir_x, ir_y_filtered, 'k', 'DisplayName', 'IR Signal'); hold on;
plot(x, y_filtered, 'r', 'DisplayName', 'Red Signal');
plot(t_ir_peaks(1:N), ir_y_filtered(locs_ir(1:N)), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'IR Peaks');
plot(t_red_peaks(1:N), y_filtered(locs_red(1:N)), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Red Peaks');
title('Bandpassed IR and Red Signals with Peak-Based Delay Detection');
xlabel('Time (s)');
ylabel('Amplitude');
legend;
grid on;

% Background noise signal
b_file = "C:\Users\keyaa\Downloads\adc_dark_LJA.csv";  % Replace with your actual filename
b_data = readmatrix(b_file);

b_x = b_data(:, 1); % better: time vector for background
b_y = b_data(:, 2); % better: background light intensity signal

Wn = [low_freq high_freq] / (fs / 2);
[b, a] = butter(order, Wn, 'Bandpass');
b_y_filtered = filtfilt(b, a, b_y); % better: remove ambient drift and noise

% Optional: Plot Background and Corrected Signal
figure;
plot(b_x, b_y_filtered, 'k');
title("Background Light Noise");
hold on;
xlabel('Time (s)');
ylabel('Amplitude (V)');
grid on;


% HEART RATE DETECTION - IR
heartbeat_times = ir_x(ir_peak_locs); % better: time points of IR heartbeats
rr_intervals = diff(heartbeat_times); % better: RR intervals in seconds
avg_rr = mean(rr_intervals); 
avg_bpm = 60 / avg_rr;
disp(['Average Heart Rate: ', num2str(avg_bpm), ' bpm']);

% HEART RATE DETECTION - Red
heartbeat_times = x(peak_locs); 
r_intervals = diff(heartbeat_times); 
avg_r = mean(r_intervals); 
avg_bpm = 60 / avg_r;
disp(['Average Red Heart Rate: ', num2str(avg_bpm), ' bpm']);


% TESTING
% Step 1: Align Red Signal to IR Using avg_delay
x_shifted = x + avg_delay;
y_red_aligned = interp1(x_shifted, y_filtered, x, 'linear', 0); % better: interpolate red onto IR time base

% Step 2: Find Peaks in Aligned Red Signal
shared_prominence = 0.01;
min_peak_dist = round(0.4 * Fs); 

[~, red_peak_locs_aligned] = findpeaks(y_red_aligned, ...
    'MinPeakProminence', shared_prominence, ...
    'MinPeakDistance', min_peak_dist);

red_peak_times_aligned = x(red_peak_locs_aligned);
rr_intervals_aligned = diff(red_peak_times_aligned);
bpm_red_aligned = 60 ./ rr_intervals_aligned;
avg_bpm_red_aligned = mean(bpm_red_aligned);
disp(['Aligned Red HR: ', num2str(avg_bpm_red_aligned), ' bpm']);

% Step 3: Compare With IR Peaks
ir_peak_times = ir_x(ir_peak_locs);
rr_ir = diff(ir_peak_times);
bpm_ir = 60 ./ rr_ir;
avg_bpm_ir = mean(bpm_ir);
disp(['IR HR: ', num2str(avg_bpm_ir), ' bpm']);

% Step 4: Plot Aligned Signals and Peaks
figure;
plot(ir_x, ir_y_filtered, 'k', 'DisplayName', 'IR Signal'); hold on;
plot(x, y_red_aligned, 'r', 'DisplayName', 'Aligned Red Signal');
plot(ir_x(ir_peak_locs), ir_y_filtered(ir_peak_locs), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'IR Peaks');
plot(x(red_peak_locs_aligned), y_red_aligned(red_peak_locs_aligned), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Red Peaks (Aligned)');
title('IR and Red Signals (Aligned) with Peaks');
xlabel('Time (s)');
ylabel('Amplitude');
legend;
grid on;

num_ir_peaks = length(ir_peak_locs);
num_red_peaks = length(peak_locs);
disp(['Number of Red Peaks (unaligned): ', num2str(num_red_peaks)]);
disp(['Number of IR Peaks:               ', num2str(num_ir_peaks)]);

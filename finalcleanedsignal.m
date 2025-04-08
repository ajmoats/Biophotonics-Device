%%% Biophotonics Checkpoint X Data Analysis:
%%% Objective: Clean and Analyze Downward Crossings of Biophotonic Signal
%%% Note: Body tremors at 3, 5, 9, 13 Hz, low frequency drift -> low-pass recommended

% Part 1: Load the Signal Data from the Oscilloscope
file = "C:\Users\keyaa\Downloads\adc_ir_AJ2_reading.csv";
data = readmatrix(file);
signal = data(:, 2);                                                        % Second .dat data is the signal amplitude
t = data(:, 1);                                                             % First .dat data is the time passed

% Part 2: Compute the Sampling Frequency and # of Samples
fs = 1 / mean(diff(t));                                                     % Sampling frequency
N = length(signal);                                                         % Total number of samples
disp(['Sampling Frequency: ', num2str(fs), ' Hz']);

% Part 3: Apply Low-pass Filter to Clean High Frequency Noise
lp_cutoff = 2;                                                              % Low-pass cutoff frequency (Hz)
lp_normalized = lp_cutoff / (fs / 2);                                       % Normalize cutoff for filter design
[b_lp, a_lp] = butter(2, lp_normalized, 'low');                             % 2nd order Butterworth low-pass
lp_signal = filtfilt(b_lp, a_lp, signal);                                   % Zero-phase filtering
disp(['Applied low-pass filter at ', num2str(lp_cutoff), ' Hz']);

% Part 4: Normalize Filtered Signal
lp_min = min(lp_signal);
lp_max = max(lp_signal);
normalized_signal = (lp_signal - lp_min) / (lp_max - lp_min);               % Normalize to [0, 1]

% Part 5: Downsample the Signal by 1000
downsampleFactor = 1000;
t_downsampled = t(1:downsampleFactor:end);
signal_downsampled = normalized_signal(1:downsampleFactor:end);

% Part 6: Smooth the Downsampled Signal (For Visual Purposes Only)
windowSize = 5;                                                             % Moving average window size
signal_smoothed = movmean(signal_downsampled, windowSize);

% Part 7: Plot the Cleaned, Downsampled, Smoothed Signal
figure;
plot(t_downsampled, signal_smoothed);
title("Low-pass Filtered, Downsampled, and Smoothed Signal");
xlabel("Time (s)");
ylabel("Normalized Amplitude");
grid on;

% Part 8: Analyze Downward Crossings at Threshold = 0.32
threshold = 0.32;
disp(['Analyzing downward crossings at threshold = ', num2str(threshold)]);

% Find where signal crosses threshold **downward**
crossings = find(signal_smoothed(1:end-1) > threshold & signal_smoothed(2:end) <= threshold);

% Print number of crossings
num_crossings = length(crossings);
disp(['Number of downward threshold crossings at ', num2str(threshold), ': ', num2str(num_crossings)]);

% Convert crossing indices to time values
crossing_times = t_downsampled(crossings);

% Apply minimum time gap to remove multiple crossings from same pulse
min_time_gap = 1.5;                                                         % Expect ~2 seconds between pulses
if ~isempty(crossing_times)
    valid_crossings = crossing_times([true; diff(crossing_times) > min_time_gap]);

    % Compute intervals between valid crossings
    if length(valid_crossings) >= 2
        intervals = diff(valid_crossings);                                  % Time between crossings (seconds)

        % Print interval and heart rate estimation
        disp('Intervals between valid downward crossings (seconds):');
        disp(intervals);

        average_interval = mean(intervals);
        estimated_hr = 60 / average_interval;
        disp(['Average interval: ', num2str(average_interval), ' seconds']);
        disp(['Estimated Heart Rate: ', num2str(estimated_hr), ' BPM']);

        % Optional: Plot crossings on signal
        figure;
        plot(t_downsampled, signal_smoothed);
        hold on;
        yline(threshold, '--r', 'Threshold 0.32');
        scatter(valid_crossings, threshold * ones(size(valid_crossings)), 'ro', 'filled');
        title("Detected Downward Crossings on Cleaned Signal");
        xlabel("Time (s)");
        ylabel("Normalized Amplitude");
        grid on;

    else
        disp('Not enough valid crossings detected for interval analysis.');
    end
else
    disp('No crossings detected.');
end

% Part 7b: Compare Full-Resolution vs. Downsampled Signal
figure;
subplot(2, 1, 1);
plot(t, normalized_signal, 'b');
title("Full-Resolution Low-Pass Filtered Signal");
xlabel("Time (s)");
ylabel("Normalized Amplitude");
grid on;

subplot(2, 1, 2);
plot(t_downsampled, signal_smoothed, 'r');
title("Downsampled and Smoothed Signal");
xlabel("Time (s)");
ylabel("Normalized Amplitude");
grid on;

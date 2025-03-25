#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import time
from gpiozero import MCP3008  # ADC
import IPython.display as display
from IPython.core.display import clear_output


ANALOG_CHANNEL = 0  # Eventually will have many
device_output = MCP3008(ANALOG_CHANNEL, max_voltage=5.0)

start_time = time.time()  # Record start time

while True:
    uptime = time.time() - start_time  # Calculate uptime
    clear_output(wait=True)  # Clears the terminal (use "cls" if on Windows)
    print(f"Uptime: {uptime:.0f} seconds")
    print(f"{device_output.value:.5f}")
    time.sleep(1)


# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')
get_ipython().run_line_magic('matplotlib', 'widget')


# In[ ]:


import time
import numpy as np
import matplotlib.pyplot as plt
from gpiozero import MCP3008

# Set up the ADC on channel 0 (scaling the reading to 0-5V)
adc = MCP3008(channel=0, max_voltage=5.0)

# Sampling parameters
sampling_interval = 0.1  # seconds per sample (roughly 20 Hz,0.05)
fs = 1 / sampling_interval  # Sampling frequency

# Create a figure with two subplots: one for time-domain and one for frequency-domain data
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Time-domain plot setup
line1, = ax1.plot([], [], lw=2)
ax1.set_ylim(0, 5)
ax1.set_xlim(0, 10)
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Voltage (V)")
ax1.set_title("Time Domain Signal")

# Frequency-domain (FFT) plot setup
line2, = ax2.plot([], [], lw=2)
ax2.set_xlabel("Frequency (Hz)")
ax2.set_ylabel("Magnitude")
ax2.set_title("Frequency Spectrum")
# Set initial x-axis up to the Nyquist frequency (half of sampling frequency)
ax2.set_xlim(0.5, fs / 2)
ax2.set_ylim(0, 0.3)  # initial limits; will adjust based on FFT magnitude

# Data buffers for time and voltage values
xdata = []
ydata = []
start_time = time.time()

try:
    while True:
        # Get current time and ADC reading
        current_time = time.time() - start_time
        voltage = adc.value * 5.0  # Scale ADC reading from 0-1 to 0-5V

        # Append new data points
        xdata.append(current_time)
        ydata.append(voltage)

        # Keep only the last 10 seconds of data for display
        while xdata and (current_time - xdata[0] > 10):
            xdata.pop(0)
            ydata.pop(0)

        # --- Update Time-Domain Plot ---
        line1.set_data(xdata, ydata)
        ax1.set_xlim(max(0, current_time - 10), current_time)

        # --- Compute and Update FFT (Frequency-Domain Plot) ---
        if len(ydata) > 1:
            y_array = np.array(ydata)
            N = len(y_array)
            # Compute one-sided FFT and corresponding frequency bins
            fft_vals = np.fft.rfft(y_array)
            fft_freq = np.fft.rfftfreq(N, d=sampling_interval)
            magnitude = np.abs(fft_vals) / N

            line2.set_data(fft_freq, magnitude)
            ax2.set_xlim(0.5, fs / 2)
            #ax2.set_ylim(0, magnitude.max() * 1.1 if magnitude.max() > 0 else 1)

        # Redraw the figure
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(sampling_interval)
except KeyboardInterrupt:
    print("Oscilloscope stopped.")


# In[ ]:





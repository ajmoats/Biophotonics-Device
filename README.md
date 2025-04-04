# Biophotonics-Device
Disease Focus: Diabetes, Hypertension, Arterial Compliance: LJ Gonzales, Omar Shehab El-Din, Alexis Moats

This GitHub supports the Spring 2025 Biophotonics Project, defined as follows:



### Introduction:
As of 2022, diabetes has risen to over 800 million cases globally, affecting 1 in 10 individuals. Individuals that have high blood pressure are twice as likely to develop diabetes. These individuals exhibit a blood pressure greater than 130/80 mmHg.

For this project we will use two wavelengths of light (∼650 and 930 nm) to measure the changes in concentration of deoxygenated hemoglobin and oxygenated hemoglobin in time. The phase difference between the two sinusoidal signals can be related to metabolism, in units of oxygen consumption per unit time. This is because the delay between the two waves tells us the speed of absorption of oxygen by cell tissue. In this project we assume that this speed of consumption of oxygen is directly related to the Metabolic Index, which indirectly reflects the consumption of glucose. We assume a linear relationship MI(t)∝d(blood glucose level)dt.

### Methods:
Data Capture: Firstly, we will measure a participant’s resting heart rate. Then, we will have participants consume a pre-workout of choice to heighten their heart rate. Over the 30 minutes, their metabolic rate will be captured by obtaining the levels of oxygenated and deoxygenated blood, extrapolating the metabolic index (MI). Then, participants will perform a workout of choice, maintaining a high heart rate, but decreasing oxygenated blood level. This will allow us to analyze the blood flow, showing overall arterial compliance.

### Structural Design: 
This cardiovascular measurement device will be a wearable ring with an internal fPCB composed of 2 sets of 2 LEDs [Figure 1]. Each of these LEDs will have 2 detectors that measure the transmission of the near infrared (NIR) scattering through the tissue. We chose a ring shape since we need to maximize surface area utilized on a thinner tissue. This in its prototype form will be connected via external wiring for now, but may utilize Bluetooth information transfer in the future. \

### Connecting Rπ to the Internet

To connect Rπ to the internet without using the HMDI-based setup, we used the bootloader to give it an initial 1) known WIFI network (mobile phone's WLAN, Author 3's MI 9T, password= [ ]). The IP address of Rπ on this network is 

### SSH-ing into Raspberry Pi and plotting data- Walkthrough
In Terminal, PASTE & ENTER:

ssh mypi@10.203.21.144

password: peanutbutter (Figure 1).

We connected (s)ecurely into the (sh)ell of Rπ through to its IP's address on the Hopkins network. This address might change if it is re-assigned (Figure 2). We now need to run the notebook locally and pass the address to the other device for visualization.

In mypi's terminal:

jupyter notebook --no-browser (without no-browser if viewing from VNC or HDMI)

In computer's terminal:

ssh -L 8888:localhost:8888 mypi@10.203.21.144 (if 8888 was one of the addresses listed by the R's jupyter kernel)

Now the local address http://localhost:8888/?token=xyz on the computer's browser should reflect that of the raspberry pi.

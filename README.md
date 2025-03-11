# Biophotonics-Device
Disease Focus: Diabetes, Hypertension, Arterial Compliance: LJ Gonzales, Omar Shehab El-Din, Alexis Moats

This GitHub supports the Spring 2025 Biophotonics Project, defined as follows:



### Introduction:
As of 2022, diabetes has risen to over 800 million cases globally, affecting 1 in 10 individuals. Individuals that have high blood pressure are twice as likely to develop diabetes. These individuals exhibit a blood pressure greater than 130/80 mmHg.

For this project we will use two wavelengths of light (∼650 and 930 nm) to measure the changes in concentration of deoxygenated hemoglobin and oxygenated hemoglobin in time. The phase difference between the two sinusoidal signals can be related to metabolism, in units of oxygen consumption per unit time. This is because the delay between the two waves tells us the speed of absorption of oxygen by cell tissue. In this project we assume that this speed of consumption of oxygen is directly related to the Metabolic Index, which indirectly reflects the consumption of glucose. We assume a linear relationship MI(t)∝d(blood glucose level)dt.

### Methods:
Data Capture: Firstly, we will measure a participant’s resting heart rate. Then, we will have participants consume a pre-workout of choice to heighten their heart rate. Over the 30 minutes, their metabolic rate will be captured by obtaining the levels of oxygenated and deoxygenated blood, extrapolating the metabolic index (MI). Then, participants will perform a workout of choice, maintaining a high heart rate, but decreasing oxygenated blood level. This will allow us to analyze the blood flow, showing overall arterial compliance.

### Structural Design: 
This cardiovascular measurement device will be a wearable ring with an internal fPCB composed of 2 sets of 2 LEDs [Figure 1]. Each of these LEDs will have 2 detectors that measure the transmission of the near infrared (NIR) scattering through the tissue. We chose a ring shape since we need to maximize surface area utilized on a thinner tissue. This in its prototype form will be connected via external wiring for now, but may utilize Bluetooth information transfer in the future. 

# Calibration Process

This directory contains scripts used to calibrate a loudspeaker array setup. The main objective of calibration is to ensure that all loudspeakers can produce a flat frequency response and are aligned to the same Sound Pressure Level (SPL) at the listening position.

## Main Calibration Script
The primary entry point is `calibration_LSRoom.m`. Running this script performs the following core steps:

1. **Microphone Calibration (Optional):** 
   You will have the option to calibrate the reference microphone using a calibrator. A 5-second recording computes an `iFactor`, which maps the digital microphone signal at 1 kHz to a standard SPL reference value.
2. **Transfer Function Measurement (`getTF.m`):**
   This script measures the Room Impulse Response (RIR) and frequency-dependent corrections for each loudspeaker by playing an exponential sweep. It outputs smoothed frequency filters to correct the response of each speaker.
3. **Level Alignment (`getLevel.m`):**
   A bandpass filtered excitation signal (Sweep, Pink Noise, or LTASS) is routed to each loudspeaker. The script measures the SPL and adjusts the loudspeaker's output level (`levelFactor`) until the targeted SPL (e.g., 70 dB) is measured within a specified tolerance.

The final calibration data is saved in two MAT-files:`current_Calibration.mat` and `currentCalibration_YYYY-MM-DD.mat`.

## Excitation Signals
When performing level alignment using LTASS (Long-Term Average Speech Spectrum), the script will prompt you to select a `.wav` file. The standard test signals (e.g., `HINT_EnglishUS_Speech_L1_S1.wav`, `HINT_Danish_Noise_Female1_Inf_new.wav`) are located in the main `src/wavFiles/` directory. The scripts are configured to point directly to this folder by default to prevent duplicate files and maintain an organized structure.

## Applying Calibration

The generated `.mat` file is loaded into the `configurationSetup` struct in the main application. 

You can apply the calibration to a given audio signal using:
- **`calibrate_vbap.m`**: Applies frequency filters and level alignments to an audio signal being panned using Vector Base Amplitude Panning (VBAP).
- **`calibrate_ambisonics.m`**: Applies frequency filters and level alignments to an Ambisonics signal before decoding to the loudspeakers.

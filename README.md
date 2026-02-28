# Iceberg: Hybrid Auralization for Minimalist Setups


**Iceberg** is a hybrid spatial audio engine designed to deliver clinical-grade virtualization using minimal hardware. By merging the localization precision of **VBAP** with the immersion of **Ambisonics**, Iceberg allows researchers to create realistic, "ecological" acoustic environments with as few as four loudspeakers.

<img width="1408" height="768" alt="Iceberg" src="https://github.com/user-attachments/assets/414c3aa9-826e-40bd-be7d-374b7da9e78c" />


<img src="https://github.com/aguirreSL/Iceberg/assets/61785053/f26b60c4-9369-4da8-a5e6-7194abb0f5cb" width="350" alt="Iceberg Auralization"><img src="https://github.com/user-attachments/assets/919a9490-44fd-4640-aeec-b41d367de4e8" width="550" alt="Iceberg Flowchart">

---

## The Problem
Individuals with hearing loss struggle in noisy, social environments. Testing their performance in these scenarios usually demands a complex laboratory. Standard virtualization methods like Ambisonics or VBAP have trade-offs. Ambisonics offers great immersion but lacks sharp localization at low orders. VBAP provides precise directionality but often lacks the "feel" or spaciousness of a real room.

## The Solution
Iceberg combines the strengths of both. It splits the Room Impulse Response (RIR) into two distinct parts:

1.  **Localization (Direct Sounde and Early Part):** Handled by **VBAP**. This ensures the listener can pinpoint exactly where a sound starts.
2.  **Immersion (Late Part):** Handled by **Ambisonics**. This provides the lush, reverberant "feel" of a real room.

The transition point isn't arbitrary. Iceberg uses **Center Time ($T_s$)** to find the ideal moment to switch from one method to the other.



---

## Key Features
* **Minimal Hardware:** Optimized for small reproduction systems (4+ speakers).
* **Dual-Listener Support:** Validated for up to two participants without losing spatial accuracy.
* **Clinically Validated:** Proven effective in speech-in-noise tasks involving both normal-hearing and hearing-aided participants.
* **High Accuracy:** Maintains localization cues within a 30° ambiguity angle in the horizontal plane.

## How it Works
The system uses MATLAB to process signals through a calibrated setup. By using $T_s$ to divide the impulse response, the early reflections that dictate spatial origin are panned precisely. The remaining energy is encoded via first-order Ambisonics to fill the room.



---

## Quick Start (For Developers)

### 1. Prerequisites
This project requires MATLAB and several specialized toolboxes. Ensure the following are in your path:
* [ITA Toolbox](https://git.rwth-aachen.de/ita/toolbox)
* [Higher-Order-Ambisonics](https://github.com/polarch/Higher-Order-Ambisonics)
* [Spherical-Harmonic-Transform](https://github.com/polarch/Spherical-Harmonic-Transform/tree/master)
* [Vector-Base-Amplitude-Panning](https://github.com/polarch/Vector-Base-Amplitude-Panning)
* [SAP Voicebox](https://github.com/ImperialCollegeLondon/sap-voicebox.git)

### 2. Implementation
* **Calibration:** Navigate to `/Calibration` and run `calibration_LSRoom.m`. Update this with your specific loudspeaker coordinates.
* **Generate Audio:** Run `iceberg_example.m` to see how to take a dry signal and an RIR to create a spatialized output.

---

## Thesis & Research
This method was developed as part of a deep dive into listening effort and virtualization. 
[**Download the full Thesis (PDF)**](https://eprints.nottingham.ac.uk/72207/1/SergioAguirre.pdf)

### Key Chapters
* **Chapter 3: Binaural Cue Distortions** – A comparison of VBAP and Ambisonics through a calibrated setup, examining spatial distortions and the impact on a second listener.
* **Chapter 4: Behavioral Study** – An investigation into how signal-to-noise ratio (SNR) and reverberation impact listening effort, using EEG and subjective questionnaires.
* **Chapter 5: The Iceberg Method** – The formal proposal and evaluation of the hybrid method using objective parameters and hearing aid verification.

### Core Findings
* **Center Time ($T_s$):** Successfully identifies the transition point between early and late reflections to split impulse responses.
* **Immersion vs. Accuracy:** The hybrid method matches the localization accuracy of VBAP while maintaining the sense of immersion typically only found in Ambisonics.
* **Clinical Viability:** The setup provides reliable binaural cues within a 30° ambiguity angle, making it suitable for audiological tests in smaller clinic rooms.

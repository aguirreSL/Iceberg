# Project Status and Todo List

## Current Status

The "Iceberg" project implements a hybrid auralization method combining VBAP (for direct sound and early reflections) and Ambisonics (for late reverberation).

### Core Components
- **Main Logic:** `src/iceberg.m` orchestrates the signal processing, splitting the impulse response and applying VBAP/Ambisonics accordingly.
- **Core Processing:** `src/iceberg_core.m` handles the splitting of the Impulse Response (IR) into early and late components based on Center Time (T20).
- **Dependencies:** The project relies on several external toolboxes:
  - ITA-Toolbox
  - VBAP (polarch)
  - HOA (polarch)
  - sap-voicebox
- **Setup:** `src/setToolboxes.m` is available to automatically check and install these dependencies.
- **Calibration:** Calibration scripts are located in `src/Calibration/`, but the integration in `iceberg.m` is currently commented out.
- **Example:** `src/iceberg_example.m` provides a usage example, demonstrating how to set up and run the auralization.

## Todo List

- [ ] **Dependency Verification**: Verify that `src/setToolboxes.m` correctly installs and configures all dependencies on a clean environment.
- [ ] **Calibration Integration**: Uncomment and verify the calibration calls in `src/iceberg.m` (`calibrate_vbap` and `calibrate_ambisonics`). Ensure the calibration logic is robust.
- [ ] **Testing**:
    - Add unit tests for `iceberg_core.m` to verify the splitting logic (e.g., check energy conservation, correct split point).
    - Create a test suite that mocks the ITA objects to verify logic without needing the full toolbox installed (if possible), or ensures the toolbox is present for tests.
- [ ] **Documentation**:
    - Improve comments in `src/iceberg.m` regarding the structure of `configSetup`.
    - Document the expected input formats for `signal` and `IR`.
- [ ] **Example Script**: Ensure `src/iceberg_example.m` runs out-of-the-box (after dependency installation).
- [ ] **Refactoring**:
    - Check if `iceberg_signal.dimensions` usage in `iceberg.m` is correct (it's used as an index limit).
    - Review error handling in `iceberg_core.m` (currently uses `try-catch` for T20 calculation).

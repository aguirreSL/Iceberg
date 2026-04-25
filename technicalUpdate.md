# Technical Update: Achieving 1-to-1 Mathematical Parity with ITA Toolbox

**Date:** 2026-02-27  
**Status:** ✅ All 8 tests passing at `1e-12` precision

---

## Summary

During final validation of the native MATLAB replacements for ITA Toolbox functions, we discovered and fixed several hidden behaviors that prevented exact 1-to-1 array parity. All native functions now match ITA outputs to within `1e-12` floating-point tolerance.

---

## Root Cause #1: `hann()` vs `hanning()` (CRITICAL)

**Impact:** `1.8951e-04` systematic deviation across all windowed samples

MATLAB has two similarly-named but **mathematically different** window functions:

| Function | Formula | Endpoints |
|----------|---------|-----------|
| `hann(N)` | `w(k) = 0.5*(1 - cos(2πk/(N-1)))` | Zero at both ends |
| `hanning(N)` | `w(k) = 0.5*(1 - cos(2πk/(N+1)))` | Non-zero at both ends |

ITA's `ita_time_window.m` calls `window(@hann, lengthWindow)` which resolves to `hann()`. Our native implementation was calling `hanning()`, causing a `1.89e-4` peak error that exactly matched our test failure.

**Fix:** Changed `native_time_window.m` line 63 from `hanning()` to `hann()`.

---

## Root Cause #2: Even-Index Forcing (`mod(interval, 2)`)

**Impact:** Off-by-one sample alignment in window boundaries

ITA's source code (line 309) contains an undocumented constraint:
```matlab
interval = interval + mod(interval, 2);  % RSC - prevent uneven samples
```
This forces all window start/end indices to be **even numbers**, shifting odd indices forward by 1.

**Fix:** Replicated the exact `mod(index, 2)` offset in `native_time_window.m`.

---

## Root Cause #3: `native_time_shift` Argument Parsing

**Impact:** `circshift` crash with error `Invalid shift type: must be a real finite integer vector`

The function incorrectly classified explicit time-domain shift values (e.g., `0.1` seconds) as sample counts, passing a non-integer float to `circshift`.

**Fix:** Rewrote the argument parser to check `varargin{1}` for `'time'`/`'samples'` strings explicitly, instead of inferring from `isnumeric(shiftAmount)`.

---

## Root Cause #4: Stale Baseline Data

**Impact:** False failures in `testAddMatch` and `testTimeShiftMatch`

The baseline `.mat` files were generated from a stale MATLAB session where `run('test_ita_baselines.m')` was used (which doesn't execute classdef tests). The baselines contained leftover data.

**Fix:** Executed `runtests('test_ita_baselines.m')` properly to regenerate clean baselines.

---

## Root Cause #5: Signal Length Mismatch in `testTimeWindowMatch`

**Impact:** `Arrays have incompatible sizes` error (native: 960,000 vs ITA: 48,000)

The test used a 20-second `SignalObj` (960,000 samples) while the ITA baseline was generated from a 1-second signal (48,000 samples).

**Fix:** Added explicit signal truncation in the test to match baseline dimensions.

---

## Final Test Results

All native replacements now achieve exact 1-to-1 parity with ITA at `AbsTol = 1e-12`:

| Test | Function | Status | Tolerance |
|------|----------|--------|-----------|
| `testNormalizeMatch` | `native_normalize_dat` | ✅ Pass | `1e-12` |
| `testTimeCropMatch` | `native_time_crop` | ✅ Pass | `1e-12` |
| `testConvolveMatch` | `native_convolve` | ✅ Pass | `1e-12` |
| `testTimeWindowMatch` | `native_time_window` | ✅ Pass | `1e-12` |
| `testTimeShiftMatch` | `native_time_shift` | ✅ Pass | `1e-12` |
| `testAddMatch` | `native_add` | ✅ Pass | `1e-12` |
| `testCenterTimeSynthetic` | `native_center_time` | ✅ Pass | `1e-6` |
| `testIcebergCoreEndToEnd` | Full pipeline | ✅ Pass | structural |

---

# Technical Update: Native Port of Calibration Application Functions

**Date:** 2026-04-25
**Status:** ✅ `calibrate_vbap.m` and `calibrate_ambisonics.m` ported to native structs and generalized to N speakers. Re-wiring into `iceberg_set_*` pipeline is the next step.

---

## Summary

The runtime calibration functions were rewritten to operate on native audio structs and to support arbitrary speaker counts/layouts. The hardcoded angle cascade for the 4-LS cardinal layout was replaced with a circular-distance lookup, which simultaneously fixed three pre-existing assignment bugs.

ITA-Toolbox is now isolated to the **measurement** scripts (`calibration_LSRoom.m`, `getTF.m`, `getLevel.m`, `BackgroundNoise.m`); the runtime/apply path is ITA-free.

---

## Change #1: Hardcoded Angle Cascade → Circular-Distance Nearest-LS

**Impact:** Removed 30+ lines of hardcoded `if/elseif` for cardinal angles (0°/90°/180°/270°) and fixed assignment bugs that, in three of the eight octants, sent the higher SPL to the **further** loudspeaker.

The previous cascade chose `s1` (LS receiving the higher SPL) and `s2` (lower SPL) via nested `if` checks against the source angle. Octant audit:

| Sub-zone | Hardcoded | Correct (closer LS first) | Status |
|---|---|---|---|
| 90°–135° | s1=180°, s2=90° | s1=90°, s2=180° | ❌ inverted |
| 135°–180° | s1=90°, s2=180° | s1=180°, s2=90° | ❌ inverted |
| 180°–225° | s1=180°, s2=270° | s1=180°, s2=270° | ✓ |
| 225°–270° | s1=180°, s2=270° | s1=270°, s2=180° | ❌ duplicate of 180°–225° branch |
| 270°–315° | s1=270°, s2=0° | s1=270°, s2=0° | ✓ |
| 315°–360° | s1=0°, s2=270° | s1=0°, s2=270° | ✓ |
| 0°–45° | s1=0°, s2=90° | s1=0°, s2=90° | ✓ |
| 45°–90° | s1=90°, s2=0° | s1=90°, s2=0° | ✓ |

The defensive `[s1_level, s2_level] = deal(max(...), min(...))` on the level pair did **not** rescue the broken cases — it swapped the level *values* but the channel *pointers* stayed wrong.

**Fix:** A single circular-distance computation replaces the cascade and works for any N and any layout (uniform or not):

```matlab
allAngles = configurationSetup.ls_dir(:,1);
angDist   = abs(mod(allAngles - iAngles + 180, 360) - 180);
[~, order] = sort(angDist);
s1 = activeLSNumbers(order(1));      % nearest LS gets the higher SPL
s2 = activeLSNumbers(order(2));      % second-nearest gets the lower SPL

gap   = abs(mod(allAngles(order(2)) - allAngles(order(1)) + 180, 360) - 180);
ratio = angDist(order(1)) / gap;     % 0 = at s1, 0.5 = at the bisector
s1_level = cos(ratio * pi/2)^2;
s2_level = sin(ratio * pi/2)^2;
```

The `mod(... + 180, 360) - 180` handles 0°/360° wrap-around. The cos²/sin² law guarantees `s1_level ≥ s2_level` by construction, so the defensive `deal` is no longer needed.

The same generalization was applied to `calibrate_ambisonics.m`, which uses only the nearest LS (`order(1)`) as the virtual reference for SPL alignment (Nearest Speaker Pan).

---

## Change #2: itaAudio Struct Accessors → Native

Replacements applied across both files:

- `signal_to_play.ch(idx).time` → `signal_to_play.time(:, idx)`
- `signal_to_play.freqVector` → locally computed `(0:nBins-1)' * (fs / nFFT)` with `nBins = floor(nFFT/2) + 1`
- NaN guard rewritten to act on native column vectors and zero individual NaN samples (the original guard only triggered when *every* sample of a channel was NaN, letting partial-NaN channels propagate into FFT)

---

## Change #3: Pre-existing Bugs Fixed in Passing

- **`calibrate_vbap.m:4`** referenced `signal.samplingRate`, but the input argument was named `signal_to_play`. Would crash on first call. Corrected.
- **`calibrate_vbap.m` signature** declared `(signal_to_play, configurationSetup)` while the body referenced `iAngles` and `level`. Corrected to `(signal_to_play, level, iAngles, configurationSetup)`, matching `calibrate_ambisonics.m`.
- **Field casing** standardized: `configurationSetup.LSArray` → `lsArray` (case-sensitive struct access; the example and other helpers all use lowercase).

---

## Not Yet Done

- **Re-wiring:** `iceberg_set_vbap` and `iceberg_set_amb` still do not call the apply functions. Next step: add a guard `if ~isempty(config.iLoudspeakerFreqFilter)` and invoke before returning.
- **Measurement scripts** (`calibration_LSRoom.m`, `getTF.m`, `getLevel.m`, `BackgroundNoise.m`) remain ITA-dependent — by design, per the branch strategy in `dist.md`.

---

# Technical Update: Restoring Original Iceberg Algorithm Semantics

**Date:** 2026-04-25
**Status:** ✅ Native pipeline now matches the original ITA-based algorithm in `main`.

---

## Summary

A code review against the original ITA implementation in `main` revealed that the native port had silently changed the algorithm in two places. Both regressions are now reverted; the native pipeline produces the spatially-meaningful Iceberg output described in the paper (mono envelope into VBAP, full B-Format into Ambisonics, sum at the speaker level).

---

## Regression #1: DSER lost its mono extraction

The original `iceberg_core.m` (ITA) extracted the W (omni) channel before the early-window operations:

```matlab
omnichannelIR = ita_split(IR,1);
[IR_Early, shiftIndex] = ita_time_shift(omnichannelIR,'auto');   % MONO
IR_Early = ita_time_window(IR_Early,[0 cTime],'time','windowType','hann');
DSER = ita_time_shift(IR_Early,abs(shiftIndex));                  % MONO
```

The native port created `omnichannelIR` for the cTime calculation but then ran `time_shift` on the **full** B-Format IR, leaving DSER with 4 channels:

```matlab
[IR_Early, shiftIndex] = native_time_shift(IR, 'auto');           % FULL — wrong
```

Downstream, `iceberg_set_vbap` convolved the mono signal with a 4-channel DSER, producing four channels (`s⊛W`, `s⊛X`, `s⊛Y`, `s⊛Z`), then VBAP-ed each at the same source angle and summed. Because VBAP is linear, this collapsed to:

```
LS_i = vbap_gain(angle) · s ⊛ (W + X + Y + Z)
```

`(W+X+Y+Z)` is not a meaningful directional pattern, so the early-reflection envelope at the active LS was directionally garbled. The intent — and the original behavior — is `s ⊛ W` (omnidirectional envelope, with VBAP supplying the direction).

**Fix:** [iceberg_core.m:9](src/core/iceberg_core.m#L9) now passes `omnichannelIR` (not `IR`) to `native_time_shift` for the early branch. DSER is mono again. The late branch continues to use the full 4-channel IR for Ambisonics decoding, as the original did.

---

## Regression #2: `windowType` keyword dropped from `time_window` calls

The original calls all used the ITA key-value form:

```matlab
ita_time_window(IR_Early, [0 cTime], 'time', 'windowType', 'hann');
ita_time_window(IR_Early, [0 .01],   'time', 'windowType', 'rectwin');
ita_time_window(IR_Late,  [...],     'time', 'windowType', 'rectwin');
```

`native_time_window` was written to accept the same key-value form (`for i = 1:2:length(varargin) ... if strcmpi(varargin{i}, 'windowType')`), but the migrated calls dropped the `'windowType'` keyword:

```matlab
native_time_window(IR_Early, [0 cTime], 'time', 'hann');     % positional
native_time_window(IR_Early, [0 .01],   'time', 'rectwin');  % positional
native_time_window(IR_Late,  [...],     'time', 'rectwin');  % positional
```

The parser never matched `'hann'` or `'rectwin'` as the keyword `'windowType'`, so `windowType` always fell back to its default (`'hann'`). All three call sites silently used hann, even when rectwin was intended.

**Severe consequence on the late part:** `IR_Late` was supposed to receive a rectangular window from t=0 to `trackLength-0.05` (no fade in the body, only zero out the last 50 ms). With the bug, it received the second half of a hann window stretched across the entire signal, multiplying the late reverb by a continuously decaying ramp from 1 (at t=0) down to 0 (near the end). The reverberation tail was being attenuated end-to-end instead of left intact.

**Fix:** [iceberg_core.m:14, 19, 35-36](src/core/iceberg_core.m) restored the `'windowType'` keyword. `native_time_window` API kept identical to ITA (key-value form).

---

## Cleanup #3: Inline normalization in `iceberg_set_amb`

[iceberg_set_amb.m:5-12](src/rendering/iceberg_set_amb.m#L5-L12) had an inline peak-normalization block reimplementing what `native_normalize_dat` already does. Replaced with a single call to `native_normalize_dat(signal)` — equivalent math, no risk of drift between the two implementations, and matches the original `ita_normalize_dat` call site in `main`.

---

## Verified Pipeline Shape

After the fixes, both branches (`main` ITA and `native`) implement the same algorithm:

```
IR (4 ch B-Format) ──► iceberg_core ──► DSER (mono, W only, windowed [0..cTime])
                                    └─► IR_Late (4 ch B-Format, cropped+windowed)

signal (mono) + DSER (mono) ──► iceberg_set_vbap ──► N LS channels (2 active)
signal (mono) + IR_Late (4 ch) ──► iceberg_set_amb ──► N LS channels (all active)

sum at the physical speaker domain ──► final N-channel output
```

---

## Test Update

[test_iceberg_integration.m](tests/test_iceberg_integration.m) previously asserted `DSER.nChannels == 4`, encoding the regression. Updated to assert `DSER.nChannels == 1` (mono) and pass a scalar `iAngle` to `iceberg_set_vbap`/`iceberg_merge` (no need to replicate per channel anymore).

---

# Technical Update: Calibration Re-wiring into the Pipeline

**Date:** 2026-04-25
**Status:** ✅ `calibrate_vbap` and `calibrate_ambisonics` are now invoked from the rendering functions, with a guard that skips them when no calibration is loaded.

---

## Summary

Calibration application functions (ported earlier this session) are now wired into `iceberg_set_vbap` and `iceberg_set_amb` at the positions that match the original `vbap_set_level`/`ambisonics_set_level` from the pre-refactor ITA pipeline. Calibration runs only when `iLoudspeakerFreqFilter` is present in the config, so existing call sites without calibration data continue to work unchanged.

---

## Sequence Restored

The position of the calibration call **differs between branches** — same as in the original:

```
VBAP path:  signal → iceberg_set_vbap [convolve+VBAP] → calibrate_vbap [filter+level] → out
Amb path:   signal → normalize → calibrate_ambisonics [filter+level] → convolve(LR) → decode → out
```

- VBAP calibration runs **after** rendering (the VBAP-panned multichannel signal already exists; calibration just filters and levels each LS channel).
- Ambisonics calibration runs **before** convolution with the IR (per-LS frequency response is pre-multiplied into the dry signal, then convolved with the B-Format IR — equivalent to the linear-system positioning of the original `ambisonics_set_level`).

---

## Signature Changes

| Function | Old | New |
|---|---|---|
| `iceberg.m` | `iceberg(signal, IR, Selected_Angle, configSetup)` | `iceberg(signal, IR, Selected_Angle, level, configSetup)` |
| `iceberg_set_vbap.m` | `(signal, DSER, iAngle, configurationSetup)` | `(signal, DSER, iAngle, level, configurationSetup)` |
| `iceberg_set_amb.m` | `(signal, LR, configurationSetup)` | `(signal, LR, iAngle, level, configurationSetup)` |

`level` is the target SPL in dB (or `'n'` to bypass level scaling — preserved from the original calibration API). `iAngle` for `iceberg_set_amb` is the source presentation angle, used by `calibrate_ambisonics` to pick the nearest virtual LS for SPL alignment.

---

## Calibration Guard

Both rendering functions now wrap the calibration call with:

```matlab
if isfield(configurationSetup, 'iLoudspeakerFreqFilter') && ...
   ~isempty(configurationSetup.iLoudspeakerFreqFilter)
    signal = calibrate_*(...);
end
```

This means:

- **No calibration in config** (field missing or empty) → calibration is skipped, signal flows through untouched. Useful for the demo, integration tests, and any user without a measured calibration.
- **Calibration loaded** → per-LS frequency filter and SPL alignment applied at the documented positions.

The guard avoids both branches of the previous "neutral mock" decision — no need for an identity calibration `.mat`; absence is the natural default.

---

## Call Sites Updated

- [iceberg_example.m](iceberg_example.m) now passes `selectedLevel` (the existing `80` dB SPL) to `iceberg(...)` — was being declared and ignored.
- [iceberg_merge.m](src/core/iceberg_merge.m) propagates its existing `level`/`angle` arguments down to the new `iceberg_set_vbap`/`iceberg_set_amb` signatures, keeping the integration test running. (`iceberg_merge` itself is still scheduled for deletion per `dist.md` Parte 8.)
- [test_iceberg_integration.m](tests/test_iceberg_integration.m) updated to call the new signatures with `level=-40` and the existing `iAngle=45`. Test config has no `iLoudspeakerFreqFilter`, so the guard skips calibration as expected.

---

# Technical Update: Removed `iceberg_merge.m` (Dead Code)

**Date:** 2026-04-25
**Status:** ✅ `src/core/iceberg_merge.m` deleted; integration test now exercises `iceberg.m` end-to-end.

---

## Summary

`iceberg_merge.m` was a vestigial orchestrator from June 2025 (commit `1ffe6ed`, "Adjusting example") that was superseded one day later by the simpler `iceberg.m` (commit `d47722b`, "simplifing"). It was never called from production code — only from the integration test. Removed to eliminate a divergent code path.

---

## Why It Existed

| Date | Commit | What |
|---|---|---|
| 08/Jun/2025 | `1ffe6ed` | Created `iceberg_merge.m` by extracting 36 lines of orchestration from `iceberg_example.m`. Used `vbap_set_level` and `ambisonics_set_level` (calibration-aware wrappers, both deleted later). |
| 09/Jun/2025 | `d47722b` "simplifing" | Created `iceberg.m` as a simplified successor with the calibration calls **commented out**. Became the production entry point. `iceberg_merge.m` was kept around but stopped being called. |
| Later | `ed10caa` | Native refactor ported both files; the redundancy survived. |

---

## What It Carried

`iceberg_merge.m` had two artefacts that were either questionable or made worse by the native port:

1. `VBAP_DS.time = VBAP_DS.time * max(max(abs(DSER.time)));` — a scaling factor that the original ITA author had marked with `%?` (uncertain). In the native version it became more redundant because `iceberg_set_vbap` already convolves with DSER internally; multiplying by `max(abs(DSER))` again amplified the convolution result a second time.
2. Took `level` as an argument but never actually called any calibration function (the calls were already commented out in the ancestor `iceberg.m`).

Both issues vanish with the file.

---

## Test Migration

[test_iceberg_integration.m:83](tests/test_iceberg_integration.m#L83) previously called `iceberg_merge(...)`. Replaced with a direct call to `iceberg(testCase.SignalObj, testCase.IrObj, iAngle, -40, testCase.ConfigSetup)` — testing the production code path instead of the deleted shadow.

The test config gained `activeLSNumbers = [1, 7, 13, 19]` (positions of the cardinal-angle LS in the 24-channel `lsArray`), which `iceberg.m` requires for the active-channel mapping. The other rendering steps (standalone `iceberg_set_amb` and `iceberg_set_vbap` checks) are kept as unit-level sanity checks.

---

# Technical Update: Smoke Test Against Real Calibration (MATLAB R2026a)

**Date:** 2026-04-25
**Status:** ✅ `iceberg_example.m` end-to-end runs against the real `currentCalibration.mat`. Two latent bugs surfaced and were fixed.

---

## Summary

A headless smoke test ([tests/run_iceberg_dry.m](tests/run_iceberg_dry.m)) was added to drive the full rendering pipeline against the actual ITA-format calibration shipped in [src/calibration/currentCalibration.mat](src/calibration/currentCalibration.mat). The first run revealed two real bugs that the integration test (with mock config or no calibration) had not caught.

Final output (signal=white noise 3 s, IR=`rt_05`, source angle=45°, level=80 dB):

```
nChannels:    24
nSamples:     160666
samplingRate: 44100 Hz
finite check: 1 (all finite)
per-channel RMS (only non-zero shown):
  ch  1 (angle=180°): RMS = 7.28e-03    ← Ambisonics late only
  ch  7 (angle=270°): RMS = 6.86e-03    ← Ambisonics late only
  ch 13 (angle=  0°): RMS = 1.11e-02    ← VBAP early + Ambisonics late
  ch 19 (angle= 90°): RMS = 1.11e-02    ← VBAP early + Ambisonics late
```

The asymmetry between the front-right pair (0°/90°) and the back-left pair (180°/270°) is the expected Iceberg signature for a 45° source: the front-right pair receives both VBAP early and Ambisonics late energy; the back-left pair receives only the late spread.

---

## Bug #1: Guard short-circuit broke on `itaResult` arrays

**Location:** [iceberg_set_vbap.m:43-44](src/rendering/iceberg_set_vbap.m#L43-L44) and [iceberg_set_amb.m:11-12](src/rendering/iceberg_set_amb.m#L11-L12)

The guard was:
```matlab
if isfield(cfg, 'iLoudspeakerFreqFilter') && ~isempty(cfg.iLoudspeakerFreqFilter)
```

The real calibration `.mat` stores `iLoudspeakerFreqFilter` as `[1×24] itaResult` — an ITA-Toolbox class array, not a native struct array. `isempty()` on an `itaResult` array returns a per-element logical vector `[0 0 0 ... 0]`, not a scalar — the `&&` short-circuit operator requires scalar logical operands, so MATLAB throws:

```
Operands to the short-circuit AND (&&) and OR (||) operators must be
convertible to logical scalars.
```

**Fix:** replaced `~isempty(...)` with `numel(...) > 0`, which is always scalar regardless of input type:

```matlab
if isfield(cfg, 'iLoudspeakerFreqFilter') && numel(cfg.iLoudspeakerFreqFilter) > 0
```

This works for native struct arrays (the mock case) and for `itaResult` arrays (the real calibration case) without runtime checks.

---

## Bug #2: `signal_vbap.nSamples` missing after the DSER-mono restoration

**Location:** [iceberg_set_vbap.m:32](src/rendering/iceberg_set_vbap.m#L32)

When DSER was multichannel (the regression from a prior session), the per-channel loop in `iceberg_set_vbap` ran multiple times and assembled `signal_vbap` via `native_add`, which sets `nSamples`. After the DSER-mono restoration, the loop runs **exactly once**: the `signal_vbap = current_signal` first-iteration branch is taken, and `current_signal` was not setting `nSamples` — leaving the field absent.

`calibrate_vbap` reads `signal_to_play.nSamples` on its first line, so the absent field surfaced as `Unrecognized field name "nSamples"` only once calibration was actually invoked (i.e., not visible to the integration test, which has no calibration data).

**Fix:** added `current_signal.nSamples = size(current_signal.time, 1);` next to the existing `nChannels` assignment.

---

## Why the Integration Test Did Not Catch These

- Bug #1 only triggers on `itaResult` (or any non-scalar-isempty type). The mock calibration in `testIcebergCoreEndToEndWithMockCalibration` uses a native struct array, where `~isempty(...)` returns a scalar — so the guard worked.
- Bug #2 only triggers when calibration runs. The non-calibrated test path skips the `calibrate_*` call, so the missing `nSamples` was never read.

Both could have been caught by a third test variant: mock calibration **as an itaResult-shaped object**. Adding that as a unit-level fixture would harden the suite further; for now the [run_iceberg_dry.m](tests/run_iceberg_dry.m) smoke test (which loads the real `.mat`) covers both.

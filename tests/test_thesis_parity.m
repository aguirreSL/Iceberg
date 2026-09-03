classdef test_thesis_parity < matlab.unittest.TestCase
    % TEST_THESIS_PARITY Pins the 2026 native chain to the thesis-era ITA
    % chain. These tests exist so the semantic drift found in the review
    % (onset-referenced centre time, Hermitian EQ mirroring, max(DSER)
    % factor, ISO 3382 onset shift) cannot silently come back.

    properties
        SampleRate = 44100
        ProjectRoot
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.ProjectRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(genpath(fullfile(testCase.ProjectRoot, 'src')));
            addpath(genpath(fullfile(testCase.ProjectRoot, 'Toolboxes')));
        end
    end

    methods(Test)
        function testCenterTimeOnsetInvariance(testCase)
            % Onset-referenced Ts must not change when the IR is delayed.
            % (The pre-fix absolute-time version shifted Ts by the delay.)
            fs = testCase.SampleRate;
            n = round(0.5*fs);
            ir = zeros(n, 1);
            ir(100) = 1;
            decay = exp(-(0:n-200).' / (0.05*fs)) * 0.1;
            ir(200:end) = ir(200:end) + decay;

            s1 = struct('time', ir, 'samplingRate', fs, 'nSamples', n, 'nChannels', 1);
            s2 = struct('time', [zeros(1000,1); ir], 'samplingRate', fs, ...
                        'nSamples', n+1000, 'nChannels', 1);
            ts1 = native_center_time(s1);
            ts2 = native_center_time(s2);
            testCase.verifyEqual(ts2, ts1, 'AbsTol', 1e-12, ...
                'Ts must be invariant to leading delay (onset-referenced).');
            testCase.verifyLessThan(ts1, 0.05, ...
                'Ts of a 50 ms-tail IR must be well under the old absolute-time value.');
        end

        function testTimeShiftSamplesMode(testCase)
            % Numeric shifts in samples must not be multiplied by fs.
            fs = 1000; N = 30000;
            x.time = zeros(N,1); x.time(40) = 1;
            x.samplingRate = fs; x.nSamples = N; x.nChannels = 1;

            [y, amt] = native_time_shift(x, 'auto');
            testCase.verifyEqual(amt, -38, ...
                'auto shift uses the ISO 3382 onset (-start+1: onset 39, impulse at 40).');
            [~, pk] = max(abs(y.time(:)));
            testCase.verifyEqual(pk, 2, 'onset shift moves the onset sample to 1 (impulse to 2).');

            z = native_time_shift(y, abs(amt), 'samples');
            [~, pk2] = max(abs(z.time(:)));
            testCase.verifyEqual(pk2, 40, 'explicit samples re-shift restores the original position.');
        end

        function testEQMirrorShelf(testCase)
            % The regression table of the review as a test: a synthetic
            % +6.02 dB shelf below 5 kHz must come out of calibrate_ambisonics
            % as a +6.02 dB shelf (pre-fix it measured +3.52 dB).
            fs = 44100; N = fs; t = (0:N-1)'/fs;
            cfg = localNeutralCal(fs);
            K = 20*log10(cfg.iFactor/2e-5);
            for f0 = [500 1000 3000 12000 20000]
                s.time = sin(2*pi*f0*t); s.samplingRate = fs; s.nSamples = N; ...
                s.nChannels = 1; s.trackLength = 1;
                y = calibrate_ambisonics(s, K, 45, cfg);  % level == K -> unit gain
                g = 20*log10(sqrt(mean(y.time.^2)) / 0.5);
                expected = 20*log10(interp1(cfg.iLoudspeakerFreqFilter(13).freqVector, ...
                    cfg.iLoudspeakerFreqFilter(13).freq, f0));
                testCase.verifyEqual(g, expected, 'AbsTol', 0.15, ...
                    sprintf('EQ at %d Hz', f0));
            end
        end

        function testVbapPanLawPin(testCase)
            % Pin the CURRENT in-pair level law (squared cosine) explicitly,
            % so any future change to it is a deliberate, visible edit.
            fs = 44100; rng(1);
            sig.time = randn(fs,1); sig.samplingRate = fs; sig.nSamples = fs; ...
            sig.nChannels = 1; sig.trackLength = 1;
            DSER.time = zeros(2048,1); DSER.time(1) = 1; DSER.samplingRate = fs; ...
            DSER.nSamples = 2048; DSER.nChannels = 1; DSER.trackLength = 2048/fs;
            cfg = localFlatCal(fs);   % flat EQ: the RMS expectation has no shelf
            cfg.ls_dir = [[180 270 0 90]; zeros(1,4)]';
            cfg.lsArray = [180:15:345 0:15:165]; cfg.activeLSNumbers = [1 7 13 19];

            K = 20*log10(cfg.iFactor/2e-5);
            for srcAngle = [22.5 30]
                v = iceberg_set_vbap(sig, DSER, srcAngle, 80, cfg);
                r = sqrt(mean(v.time.^2, 1));
                rho = srcAngle/90;                       % position inside 0/90 pair
                a1 = cos(rho*pi/2)^2; a2 = sin(rho*pi/2)^2;
                testCase.verifyEqual(20*log10(r(3)/r(4)), 20*log10(a1/a2), ...
                    'AbsTol', 0.1, 'in-pair ratio must follow the cos^2/sin^2 law');
                testCase.verifyEqual(r(3), a1 * 0.5 * 10^((80-K)/20), 'RelTol', 0.05, ...
                    'active-channel RMS must equal 0.5*10^((L-K)/20)*a_i');
            end
        end

        function testMaxDSERFactor(testCase)
            % Thesis chain: VBAP_DS = VBAP_DS * max(DSER). A DSER peaking at
            % 0.707 (FuMa W) must scale the rendered pair by 0.707.
            fs = 44100; rng(2);
            sig.time = randn(fs,1); sig.samplingRate = fs; sig.nSamples = fs; ...
            sig.nChannels = 1; sig.trackLength = 1;
            DSER.time = zeros(2048,1); DSER.time(1) = 0.7071; DSER.samplingRate = fs; ...
            DSER.nSamples = 2048; DSER.nChannels = 1; DSER.trackLength = 2048/fs;
            cfg = localFlatCal(fs);
            cfg.ls_dir = [[180 270 0 90]; zeros(1,4)]';
            cfg.lsArray = [180:15:345 0:15:165]; cfg.activeLSNumbers = [1 7 13 19];
            K = 20*log10(cfg.iFactor/2e-5);

            v = iceberg_set_vbap(sig, DSER, 0, 80, cfg);   % source at the 0 deg LS
            r = sqrt(mean(v.time.^2, 1));
            testCase.verifyEqual(r(3), 0.7071 * 0.5 * 10^((80-K)/20), 'RelTol', 0.05, ...
                'output must carry the max(DSER) factor');
        end

        function testLevelCharNRunsThroughPipeline(testCase)
            % 'n' must not crash the guard (regression: 'n' > 95 compared a
            % char with a number and threw). Thesis-inherited semantics are
            % pinned per branch: amb skips the level term, vbap scales to
            % -lsdBperVolt.
            fs = 44100; rng(4); N = fs;
            s.time = randn(N,1); s.samplingRate = fs; s.nSamples = N; ...
            s.nChannels = 1; s.trackLength = 1;
            s4.time = repmat(s.time, 1, 4); s4.samplingRate = fs; s4.nSamples = N; ...
            s4.nChannels = 4; s4.trackLength = 1;
            cfg = localFlatCal(fs);
            K = 20*log10(cfg.iFactor/2e-5);

            yA = calibrate_ambisonics(s, 'n', 45, cfg);
            testCase.verifyEqual(sqrt(mean(yA.time.^2)), 0.5, 'RelTol', 1e-6, ...
                "amb branch 'n': no level term (flat EQ, Gamma=1)");

            yV = calibrate_vbap(s4, 'n', 45, cfg);
            testCase.verifyEqual(sqrt(mean(yV.time(:,3).^2)), ...
                0.5 * 10^(-K/20), 'RelTol', 1e-3, ...
                "vbap branch 'n': levels stay at 0, so output scales to -lsdBperVolt (no a_i)");

            % IR must exceed 50 ms: the late window [0, T-0.05] silently
            % empties shorter IRs (known edge behaviour of the tail window).
            fs0 = fs; N0 = round(0.6*fs0);
            ir = zeros(N0, 4); ir(1475, :) = 1; ir(2000, :) = ir(2000, :) + 0.01;
            IR.time = ir; IR.samplingRate = fs0; IR.nSamples = N0; IR.nChannels = 4; ...
            IR.trackLength = N0/fs0;
            out = iceberg(s, IR, 45, 'n', cfg);   % must not throw
            testCase.verifyTrue(all(isfinite(out.time(:))), "'n' output is finite");
        end

        function testEQOddLengthSignal(testCase)
            % Odd nFFT takes the other mirroring branch; the shelf must hold.
            fs = 44100; N = fs + 1; t = (0:N-1)'/fs;
            cfg = localNeutralCal(fs);
            K = 20*log10(cfg.iFactor/2e-5);
            for f0 = [1000 8000]
                s.time = sin(2*pi*f0*t); s.samplingRate = fs; s.nSamples = N; ...
                s.nChannels = 1; s.trackLength = N/fs;
                y = calibrate_ambisonics(s, K, 45, cfg);
                g = 20*log10(sqrt(mean(y.time.^2)) / 0.5);
                expected = 20*log10(interp1(cfg.iLoudspeakerFreqFilter(13).freqVector, ...
                    cfg.iLoudspeakerFreqFilter(13).freq, f0));
                testCase.verifyEqual(g, expected, 'AbsTol', 0.2, ...
                    sprintf('odd-length EQ at %d Hz', f0));
            end
        end

        function testEQComplexFilterData(testCase)
            % The shipped .mat stores itaResult with complex freq data; the
            % mirroring must survive complex curves and keep their magnitude.
            fs = 44100; N = fs; t = (0:N-1)'/fs;
            cfg = localFlatCal(fs);
            phase = 0.3;
            for i = 1:24
                cfg.iLoudspeakerFreqFilter(i).freq = ...
                    cfg.iLoudspeakerFreqFilter(i).freq * exp(1i*phase);
            end
            K = 20*log10(cfg.iFactor/2e-5);
            s.time = sin(2*pi*1000*t); s.samplingRate = fs; s.nSamples = N; ...
            s.nChannels = 1; s.trackLength = 1;
            y = calibrate_ambisonics(s, K, 45, cfg);
            testCase.verifyTrue(all(isfinite(y.time(:))), 'complex EQ gives finite output');
            g = 20*log10(sqrt(mean(y.time.^2)) / 0.5);
            testCase.verifyEqual(g, 0, 'AbsTol', 0.2, ...
                'constant-phase complex curve keeps unit magnitude at 1 kHz');
        end

        function testStartIRPreEcho(testCase)
            % ISO 3382 anti-oscillation rule: energy excursions above
            % threshold+6 dB before the onset move the onset in front of them.
            fs = 44100; N = round(0.2*fs);
            ir = zeros(N,1); ir(50) = 10^(-12/20); ir(100) = 1;   % -12 dB pre-echo
            got = native_start_IR(struct('time', ir, 'samplingRate', fs, ...
                'nSamples', N, 'nChannels', 1));
            testCase.verifyEqual(got, 49, 'onset moves before the pre-echo burst');

            if exist('ita_start_IR', 'file') == 2
                ref = ita_start_IR(itaAudio(ir, fs, 'time'));
                testCase.verifyEqual(double(ref), got, 'must match ITA');
            end
        end

        function testStartIRBadSNR(testCase)
            % White noise is not an IR: the SNR guard returns onset 1.
            fs = 44100; rng(5);
            noise = randn(round(0.5*fs), 1);
            got = native_start_IR(struct('time', noise, 'samplingRate', fs, ...
                'nSamples', numel(noise), 'nChannels', 1));
            testCase.verifyEqual(got, 1, 'no onset shift for SNR < 20 dB');
        end

        function testTimeShiftThresholdString(testCase)
            fs = 1000; N = 20000;
            x.time = zeros(N,1); x.time(100) = 1;
            x.samplingRate = fs; x.nSamples = N; x.nChannels = 1;
            [y, amt] = native_time_shift(x, '20dB');
            testCase.verifyEqual(amt, -98, "'20dB' shift: onset 99 -> -99+1");
            [~, pk] = max(abs(y.time(:)));
            testCase.verifyEqual(pk, 2, 'impulse lands at onset+1');
        end

        function testAnechoicThroughIceberg(testCase)
            % configSetup.anechoicSpecialCase plumbs through iceberg().
            fs = 44100; N = round(0.6*fs);
            ir = zeros(N, 4); ir(1475, :) = 1; ir(2000, :) = ir(2000, :) + 0.01;
            IR.time = ir; IR.samplingRate = fs; IR.nSamples = N; IR.nChannels = 4; ...
            IR.trackLength = N/fs;
            rng(6); sig.time = randn(round(fs),1); sig.samplingRate = fs; ...
            sig.nSamples = fs; sig.nChannels = 1; sig.trackLength = 1;

            cfg = localFlatCal(fs);
            cfg.anechoicSpecialCase = true;
            out = iceberg(sig, IR, 45, 80, cfg);
            testCase.verifyEqual(out.nChannels, 24, 'full master array');
            activeRms = sqrt(mean(out.time(:, cfg.activeLSNumbers).^2, 1));
            testCase.verifyTrue(all(activeRms > 0), ...
                'both branches audible through the anechoic path');
        end

        function testAnechoicSpecialCase(testCase)
            fs = 44100; N = round(0.6*fs);
            ir = zeros(N, 4); ir(1475, :) = 1;   % bare direct peak
            ir(2000, :) = ir(2000, :) + 0.01;    % one late reflection (~12 ms)
            IR.time = ir; IR.samplingRate = fs; IR.nSamples = N; IR.nChannels = 4; ...
            IR.trackLength = N/fs;
            [D, L] = iceberg_core(IR, true);
            testCase.verifyEqual(D.time(1475), 1, 'AbsTol', 1e-12, ...
                'anechoic DSER keeps the direct peak at its absolute position');
            testCase.verifyEqual(sum(abs(D.time(1:1473))), 0, ...
                'nothing before the onset');
            testCase.verifyEqual(sum(abs(D.time(1916:end))), 0, ...
                'anechoic DSER is zero beyond onset+10 ms');
            testCase.verifyGreaterThan(sum(L.time(:).^2), 0, ...
                'anechoic late branch keeps the post-10ms remainder');
        end

        function testStartIRMatchesITA(testCase)
            testCase.assumeTrue(exist('ita_start_IR', 'file') == 2, 'ITA Toolbox not installed');
            files = {'rt_00', 'rt_05', 'rt_11'};
            for iF = 1:numel(files)
                [x, f] = audioread(fullfile(testCase.ProjectRoot, 'src', 'wavFiles', ...
                    files{iF}, 'BFormat1.Wav'));
                wi = itaAudio(x, f, 'time');
                ref = ita_start_IR(wi);
                got = native_start_IR(struct('time', x, 'samplingRate', f, ...
                    'nSamples', size(x,1), 'nChannels', size(x,2)));
                testCase.verifyEqual(got(:)', double(ref(:))', ...
                    sprintf('onset indices must match ITA (%s)', files{iF}));
            end
        end

        function testCenterTimeMatchesITA(testCase)
            testCase.assumeTrue(exist('ita_roomacoustics', 'file') == 2, 'ITA Toolbox not installed');
            files = {'rt_05', 'rt_11'};
            for iF = 1:numel(files)
                [x, f] = audioread(fullfile(testCase.ProjectRoot, 'src', 'wavFiles', ...
                    files{iF}, 'BFormat1.Wav'));
                wi = itaAudio(x(:,1), f, 'time');
                c = ita_roomacoustics(wi, 'Center_Time', 'broadbandAnalysis', 1);
                ctRef = double(c.Center_Time.freq);
                if isnan(ctRef)
                    c = ita_roomacoustics(ita_normalize_dat(wi), 'Center_Time', ...
                        'broadbandAnalysis', 1, 'edcMethod', 'noCut');
                    ctRef = double(c.Center_Time.freq);
                end
                ctGot = native_center_time(struct('time', x(:,1), 'samplingRate', f, ...
                    'nSamples', size(x,1), 'nChannels', 1));
                testCase.verifyEqual(ctGot, ctRef, 'AbsTol', 5e-4, ...
                    sprintf('Ts must match ITA (%s)', files{iF}));
            end
        end

        function testEQMatchesItaMultiplySpk(testCase)
            % Sharpest EQ check: same input, same curve, native mirrored
            % multiplication vs ita_multiply_spk, sample by sample.
            testCase.assumeTrue(exist('ita_multiply_spk', 'file') == 2, 'ITA Toolbox not installed');
            fs = 44100; rng(3); N = 2*fs;
            x = randn(N,1); xn = x / (2*sqrt(mean(x.^2)));
            fv = linspace(0, fs/2, 1024)';
            H = ones(1024,1); H(fv < 5000) = 1.25; H(fv > 15000) = 0.8;

            cfg = struct();
            cfg.ls_dir = [[180 270 0 90]; zeros(1,4)]';
            cfg.lsArray = [180:15:345 0:15:165]; cfg.activeLSNumbers = [1 7 13 19];
            cfg.iFactor = 1; cfg.newLevelFactor = ones(1,24);
            for i = 1:24
                cfg.iLoudspeakerFreqFilter(i).freqVector = fv;
                cfg.iLoudspeakerFreqFilter(i).freq = H;
            end
            s.time = xn; s.samplingRate = fs; s.nSamples = N; s.nChannels = 1; s.trackLength = N/fs;
            K = 20*log10(cfg.iFactor/2e-5);
            yNat = calibrate_ambisonics(s, K, 45, cfg);   % unit level gain, Gamma=1

            Interp = zeros(numel(itaAudio(xn, fs, 'time').freqVector), 19);
            Interp(:,13) = pchip(fv, H, localFreqVector(xn, fs));
            filt = itaAudio(Interp, fs, 'freq');
            yIta = ita_multiply_spk(itaAudio(xn, fs, 'time'), filt.ch(13));

            testCase.verifyEqual(yNat.time(:), yIta.time(:), 'RelTol', 1e-10, ...
                'native EQ multiplication must match ita_multiply_spk');
        end

        function testPairSelectionModes(testCase)
            % At 100 deg the 2022 cascade sent the loud coefficient to the
            % 180 deg loudspeaker (inversion); the corrected selection sends
            % it to 90 deg. Default must be the thesis behaviour (bit
            % parity, error included); 'nearest' keeps the 25 Apr 2026 fix.
            fs = 44100; rng(12);
            sig.time = randn(fs,1); sig.samplingRate = fs; sig.nSamples = fs; ...
            sig.nChannels = 1; sig.trackLength = 1;
            DSER.time = zeros(2048,1); DSER.time(1) = 1; DSER.samplingRate = fs; ...
            DSER.nSamples = 2048; DSER.nChannels = 1; DSER.trackLength = 2048/fs;
            cfg = localFlatCal(fs);
            cfg.ls_dir = [[180 270 0 90]; zeros(1,4)]';
            cfg.lsArray = [180:15:345 0:15:165]; cfg.activeLSNumbers = [1 7 13 19];

            vCompat = iceberg_set_vbap(sig, DSER, 100, 80, cfg);
            rC = sqrt(mean(vCompat.time.^2, 1));
            testCase.verifyLessThan(rC(3), rC(1), ...
                'cascade2022 (default): 100 deg source is louder on the 180 deg LS (thesis error reproduced)');

            cfg.pairSelection = 'nearest';
            vNear = iceberg_set_vbap(sig, DSER, 100, 80, cfg);
            rN = sqrt(mean(vNear.time.^2, 1));
            testCase.verifyLessThan(rN(1), rN(4), ...
                'nearest mode: 100 deg source is louder on the 90 deg LS (corrected)');
        end

        function testNearestLSTieBreak(testCase)
            % At the 45/135/225/315 equidistant points the 2022 cascade
            % picked the loudspeaker CLOCKWISE from the source (45->0,
            % 135->90, 225->180, 315->270). Pin the same tie-break.
            fs = 44100; rng(8); N = fs;
            fv = linspace(0, fs/2, 512)';
            cfg = struct();
            cfg.ls_dir = [[180 270 0 90]; zeros(1,4)]';
            cfg.lsArray = [180:15:345 0:15:165]; cfg.activeLSNumbers = [1 7 13 19];
            cfg.iFactor = 1; cfg.newLevelFactor = ones(1,24);
            for i = 1:24
                cfg.iLoudspeakerFreqFilter(i).freqVector = fv;
                cfg.iLoudspeakerFreqFilter(i).freq = i * ones(512,1);  % unique gain per LS
            end
            K = 20*log10(cfg.iFactor/2e-5);
            ties = {45, 13; 135, 19; 225, 1; 315, 7};   % angle -> expected LS (physical index)
            for t = 1:size(ties,1)
                s.time = randn(N,1); s.samplingRate = fs; s.nSamples = N; ...
                s.nChannels = 1; s.trackLength = 1;
                y = calibrate_ambisonics(s, K, ties{t,1}, cfg);
                testCase.verifyEqual(sqrt(mean(y.time.^2)), 0.5 * ties{t,2}, 'RelTol', 1e-9, ...
                    sprintf('tie at %d deg must take the EQ of physical LS %d', ties{t,1}, ties{t,2}));
            end
        end

        function testEndToEndChainMatchesITA(testCase)
            % Full-pipeline golden test: iceberg() versus a thesis-chain
            % reconstruction built inline from ITA primitives (the path the
            % thesis was measured with). Uses the committed rt_05 IR and the
            % shipped calibration. Angle 45 deg: both pair-selection policies
            % (2022 cascade and nearest-LS) agree there.
            testCase.assumeTrue(exist('ita_roomacoustics', 'file') == 2, 'ITA Toolbox not installed');

            fs = 44100; Tsig = 1; lvl = 80; iAngle = 45;
            projectRoot = testCase.ProjectRoot;

            [x, f] = audioread(fullfile(projectRoot, 'src', 'wavFiles', 'rt_05', 'BFormat1.Wav'));
            Sc = load(fullfile(projectRoot, 'src', 'calibration', 'currentCalibration.mat'), ...
                'newLevelFactor', 'iFactor', 'iLoudspeakerFreqFilter');

            cfg = struct();
            cfg.lsArray = [180:15:345 0:15:165];
            cfg.ls_dir  = [[180 270 0 90]; zeros(1,4)]';
            cfg.activeLSNumbers = [1 7 13 19];
            cfg.newLevelFactor = Sc.newLevelFactor;
            cfg.iFactor = Sc.iFactor;
            for i = 1:24
                cfg.iLoudspeakerFreqFilter(i).freqVector = double(Sc.iLoudspeakerFreqFilter(i).freqVector(:));
                cfg.iLoudspeakerFreqFilter(i).freq       = Sc.iLoudspeakerFreqFilter(i).freq(:);
            end

            % ---------- reference: thesis chain with ITA primitives ----------
            rng(9); xn = randn(Tsig*fs, 1);
            sigI = ita_normalize_dat(itaAudio(xn, fs, 'time'));
            iri  = itaAudio(x, f, 'time');
            omni = ita_split(iri, 1);

            c = ita_roomacoustics(omni, 'Center_Time', 'broadbandAnalysis', 1);
            cTime = double(c.Center_Time.freq);
            if isnan(cTime)
                c = ita_roomacoustics(ita_normalize_dat(omni), 'Center_Time', ...
                    'broadbandAnalysis', 1, 'edcMethod', 'noCut');
                cTime = double(c.Center_Time.freq);
            end

            % DSER (early branch)
            [IR_Early, sIdx] = ita_time_shift(omni, 'auto');
            IR_Early = ita_time_window(IR_Early, [0 cTime], 'time', 'windowType', 'hann');
            DSERi = ita_time_shift(IR_Early, abs(sIdx), 'time');

            % amb branch: calibrate dry signal (nearest-LS EQ), convolve late, decode
            [D4, ~] = ambiDecoder(cfg.ls_dir, 'SAD', 1, 1);
            angDist = abs(mod(cfg.ls_dir(:,1) - iAngle + 180, 360) - 180);
            [~, order] = sort(angDist);
            iChannel = cfg.activeLSNumbers(order(1));
            lsdB = 20*log10(cfg.iFactor/2e-5);
            fqv = sigI.freqVector; fqv = double(reshape(fqv, [], 1));
            InterpA = zeros(numel(fqv), 19);
            InterpA(:, iChannel) = pchip(cfg.iLoudspeakerFreqFilter(iChannel).freqVector, ...
                cfg.iLoudspeakerFreqFilter(iChannel).freq, fqv);
            filtA = itaAudio(InterpA, fs, 'freq');
            xr = sigI.time(:,1);
            xr = xr / (2*sqrt(mean(xr.^2))) * 10^((lvl - lsdB)/20);
            sigAmb = ita_multiply_spk(itaAudio(xr, fs, 'time'), filtA.ch(iChannel)) * cfg.newLevelFactor(iChannel);

            [IR_LateI, s2] = ita_time_shift(iri, 'auto');
            IR_LateI = ita_time_crop(IR_LateI, [cTime 0], 'time');
            resync = iri.nSamples - IR_LateI.nSamples;
            IR_LateI.time = [zeros(resync, 4); IR_LateI.time];
            IR_LateI = ita_time_window(IR_LateI, [0 (IR_LateI.trackLength - 0.05)], 'time', 'windowType', 'rectwin');
            IR_LateI = ita_time_shift(IR_LateI, abs(s2), 'time');
            ambRef = decodeBformat(ita_convolve(sigAmb, IR_LateI).time, D4);

            % vbap branch: convolve, ring_VBAP, calibrate per channel, x max(DSER)
            convI = ita_convolve(sigI, DSERi);
            panI = local_ring_VBAP(fs, convI.time(:,1), iAngle, cfg);
            s1Db = 20*log10(10^(lvl/20) * cos((iAngle/90)*pi/2)^2);
            s2Db = 20*log10(10^(lvl/20) * sin((iAngle/90)*pi/2)^2);
            levels = zeros(1, 19);
            levels(cfg.activeLSNumbers(order(1))) = s1Db;
            levels(cfg.activeLSNumbers(order(2))) = s2Db;
            nConv = size(convI.time, 1);
            fqv2 = (0:floor(nConv/2))' * (fs / nConv);
            InterpV = zeros(numel(fqv2), 19);
            for ic = cfg.activeLSNumbers
                InterpV(:, ic) = pchip(cfg.iLoudspeakerFreqFilter(ic).freqVector, ...
                    cfg.iLoudspeakerFreqFilter(ic).freq, fqv2);
            end
            filtV = itaAudio(InterpV, fs, 'freq');
            vbRef = zeros(size(panI,1), 4);
            for ic = 1:4
                ch = panI(:, ic);
                ch = ch / (2*sqrt(mean(ch.^2)));
                ch = ch * 10^((levels(cfg.activeLSNumbers(ic)) - lsdB)/20);
                y = ita_multiply_spk(itaAudio(ch, fs, 'time'), filtV.ch(cfg.activeLSNumbers(ic)));
                vbRef(:, ic) = y.time(:,1) * cfg.newLevelFactor(cfg.activeLSNumbers(ic));
            end
            vbRef = vbRef * max(DSERi.time(:));

            % map into the master array and add
            nRef = max(size(vbRef,1), size(ambRef,1));
            ref = zeros(nRef, 24);
            for ic = 1:4
                ref(1:size(vbRef,1), cfg.activeLSNumbers(ic)) = ...
                    ref(1:size(vbRef,1), cfg.activeLSNumbers(ic)) + vbRef(:, ic);
                ref(1:size(ambRef,1), cfg.activeLSNumbers(ic)) = ...
                    ref(1:size(ambRef,1), cfg.activeLSNumbers(ic)) + ambRef(:, ic);
            end

            % ---------- native ----------
            s3.time = xn; s3.samplingRate = fs; s3.nSamples = numel(xn); ...
            s3.nChannels = 1; s3.trackLength = Tsig;
            IRs.time = x; IRs.samplingRate = f; IRs.nSamples = size(x,1); ...
            IRs.nChannels = size(x,2); IRs.trackLength = size(x,1)/f;
            outN = iceberg(s3, IRs, iAngle, lvl, cfg);

            n = min(size(ref,1), size(outN.time,1));
            for ic = cfg.activeLSNumbers
                refMax = max(abs(ref(:, ic)));
                if refMax == 0, continue; end
                d = max(abs(ref(1:n, ic) - outN.time(1:n, ic)));
                testCase.verifyLessThan(d / refMax, 1e-6, ...
                    sprintf('channel %d must match the ITA chain (rel)', ic));
            end
        end
    end
end

function pansig = local_ring_VBAP(fs, iSignal, iAngle, configurationSetup)
    % Verbatim copy of the ring_VBAP subfunction of iceberg_set_vbap
    % (the pan stage is shared verbatim between both chains).
    blocksize = fs/18.3750;
    if fs == 48000
        blocksize = fs/20;
    elseif fs == 96000
        blocksize = fs/40;
    end
    hopsize = blocksize/2;
    ls_num = length(configurationSetup.ls_dir);
    sig = iSignal;
    Lsig = length(sig);
    Nhop = ceil(Lsig/hopsize) + 2;
    padsig = [zeros(hopsize,1); sig; zeros(Nhop*hopsize - Lsig - hopsize,1)];
    pansig = zeros(size(padsig,1), ls_num);
    static = ones(length((0:(Nhop-1)-1)'*(9*360)/(Nhop-1)),1);
    iAngleCount = iAngle;
    azis = iAngleCount*static;
    eles = 0*azis;
    ls_groups = findLsPairs(configurationSetup.ls_dir(:,1));
    layoutInvMtx = invertLsMtx(configurationSetup.ls_dir(:,1), ls_groups);
    counter = 1;
    window = hanning(blocksize);
    spread = 0;
    for idx = 0:hopsize:(Nhop-2)*hopsize
        winsig = padsig(idx+(1:blocksize),1).*window;
        azi = azis(counter);
        gains = vbap([azi 0], ls_groups, layoutInvMtx, spread);
        panwinsig = winsig*gains;
        pansig(idx+(1:blocksize),:) = pansig(idx+(1:blocksize),:) + panwinsig;
        counter = counter+1;
    end
    pansig = pansig(hopsize+(1:Lsig),:);
end

function fv = localFreqVector(x, fs)
    nBins = floor(numel(x)/2) + 1;
    fv = (0:nBins-1)' * (fs / numel(x));
end

function cfg = localFlatCal(fs)
    nBins = 1024;
    fv = linspace(0, fs/2, nBins)';
    H = ones(nBins, 1);
    cfg = struct();
    cfg.ls_dir = [[180 270 0 90]; zeros(1,4)]';
    cfg.lsArray = [180:15:345 0:15:165]; cfg.activeLSNumbers = [1 7 13 19];
    cfg.iFactor = 1; cfg.newLevelFactor = ones(1,24);
    for i = 1:24
        cfg.iLoudspeakerFreqFilter(i).freqVector = fv;
        cfg.iLoudspeakerFreqFilter(i).freq = H;
    end
end

function cfg = localNeutralCal(fs)
    nBins = 1024;
    fv = linspace(0, fs/2, nBins)';
    H = ones(nBins, 1); H(fv < 5000) = 2;      % +6.02 dB below 5 kHz
    H(fv > 15000) = 0.794328234724282;          % -2 dB above 15 kHz
    cfg = struct();
    cfg.ls_dir = [[180 270 0 90]; zeros(1,4)]';
    cfg.lsArray = [180:15:345 0:15:165]; cfg.activeLSNumbers = [1 7 13 19];
    cfg.iFactor = 1; cfg.newLevelFactor = ones(1,24);
    for i = 1:24
        cfg.iLoudspeakerFreqFilter(i).freqVector = fv;
        cfg.iLoudspeakerFreqFilter(i).freq = H;
    end
end

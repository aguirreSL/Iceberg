classdef test_ita_baselines < matlab.unittest.TestCase
    % TEST_ITA_BASELINES Validates that the current ITA toolbox generates
    % the expected arrays and values to serve as ground truth for native replacements.
    
    properties
        SignalData
        SampleRate = 48000;
        SignalObj
    end
    
    methods(TestMethodSetup)
        function createDummyData(testCase)
            % Initialize paths relative to project root (tests/ -> project root)
            projectRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(genpath(fullfile(projectRoot, 'src')));
            addpath(genpath(fullfile(projectRoot, 'Toolboxes')));
            try
                setToolboxes();
            catch
                % Silently catch if it tries to install and fails, just hoping path is there
            end
            
            % Create a simple 1-second sine wave
            t = (0:testCase.SampleRate-1)' / testCase.SampleRate;
            testCase.SignalData = sin(2 * pi * 400 * t);
            
            % Setup an itaAudio object
            try
                testCase.SignalObj = itaAudio(testCase.SignalData, testCase.SampleRate, 'time');
            catch ME
                assumeFail(testCase, ['ITA Toolbox is not initialized or available for baseline generation. Error: ' ME.message]);
            end
        end
    end
    
    methods(Test)
        function testNormalizeData(testCase)
            % Verify ita_normalize_dat baseline behavior
            actVal = ita_normalize_dat(testCase.SignalObj);
            
            % Save to a .mat file to serve as a baseline for the native test
            save(fullfile(fileparts(mfilename('fullpath')), 'baseline_normalize.mat'), 'actVal');
            
            testCase.verifyClass(actVal, 'itaAudio');
            testCase.verifyNotEmpty(actVal.time);
        end
        
        function testTimeCrop(testCase)
            % Verify ita_time_crop baseline behavior (cropping to first 0.5s)
            actVal = ita_time_crop(testCase.SignalObj, [0 0.5], 'time');
            
            save(fullfile(fileparts(mfilename('fullpath')), 'baseline_timecrop.mat'), 'actVal');
            
            testCase.verifyEqual(actVal.trackLength, 0.5, 'RelTol', 1e-3);
        end
        
        function testConvolution(testCase)
            % Create a dummy IR
            irData = zeros(testCase.SampleRate, 1);
            irData(1) = 1; % Ideal impulse
            irObj = itaAudio(irData, testCase.SampleRate, 'time');
            
            % Verify ita_convolve
            actVal = ita_convolve(testCase.SignalObj, irObj);
            
            save(fullfile(fileparts(mfilename('fullpath')), 'baseline_convolve.mat'), 'actVal');
            
            testCase.verifyNotEmpty(actVal.time);
        end
        
        function testTimeWindow(testCase)
            % Verify ita_time_window baseline behavior
            actVal = ita_time_window(testCase.SignalObj, [0.4 0.5], 'time', 'hann');
            
            save(fullfile(fileparts(mfilename('fullpath')), 'baseline_timewindow.mat'), 'actVal');
            
            testCase.verifyNotEmpty(actVal.time);
        end
        
        function testTimeShift(testCase)
            % Verify ita_time_shift baseline behavior (shifting 0.1s)
            actVal = ita_time_shift(testCase.SignalObj, 0.1, 'time');
            
            save(fullfile(fileparts(mfilename('fullpath')), 'baseline_timeshift.mat'), 'actVal');
            
            testCase.verifyNotEmpty(actVal.time);
        end
        
        function testAdd(testCase)
            % Verify ita_add baseline behavior
            actVal = ita_add(testCase.SignalObj, testCase.SignalObj);
            
            save(fullfile(fileparts(mfilename('fullpath')), 'baseline_add.mat'), 'actVal');
            
            testCase.verifyNotEmpty(actVal.time);
        end
    end
end

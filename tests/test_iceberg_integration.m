classdef test_iceberg_integration < matlab.unittest.TestCase
    % TEST_ICEBERG_INTEGRATION Validates mapping of the final iceberg core integration.
    
    properties
        SignalObj
        SampleRate = 48000;
        IrData
        IrObj
        ConfigSetup
    end
    
    methods(TestMethodSetup)
        function loadGroundTruth(testCase)
            % Initialize paths and tools
            addpath(genpath(fullfile(pwd, 'src')));
            addpath(genpath(fullfile(pwd, 'Toolboxes')));
            try
                setToolboxes();
            catch
                % Silently catch
            end
            
            % Create the identical base signals used in baseline generation
            t = (0:testCase.SampleRate-1)' / testCase.SampleRate;
            signalData = sin(2 * pi * 400 * t);
            
            % Mock the structure of native audioStruct
            testCase.SignalObj = struct();
            testCase.SignalObj.time = signalData;
            testCase.SignalObj.samplingRate = testCase.SampleRate;
            testCase.SignalObj.nSamples = length(signalData);
            testCase.SignalObj.nChannels = 1;
            testCase.SignalObj.trackLength = testCase.SignalObj.nSamples / testCase.SampleRate;
            
            % Mock a simple impulse
            testCase.IrData = zeros(testCase.SampleRate, 4); % 4 channel B-Format mock
            testCase.IrData(100, :) = 1; % Ideal impulse delayed by a fraction
            
            testCase.IrObj = struct();
            testCase.IrObj.time = testCase.IrData;
            testCase.IrObj.samplingRate = testCase.SampleRate;
            testCase.IrObj.nSamples = length(testCase.IrData);
            testCase.IrObj.nChannels = 4;
            testCase.IrObj.trackLength = testCase.IrObj.nSamples / testCase.SampleRate;
            
            % Mock configuration structure
            testCase.ConfigSetup = struct();
            iceberglAngles = [180, 270, 0, 90]; 
            testCase.ConfigSetup.ls_dir = [iceberglAngles; zeros(1,4)]';
            testCase.ConfigSetup.lsArray = [180,195,210,225,240,255,270,285,300,315,330,345, ...
                            0,15,30,45,60,75,90,105,120,135,150,165];
        end
    end
    
    methods(Test)
        function testIcebergCoreEndToEnd(testCase)
            % This simulates replacing `iceberg.m` with the native logic
            
            % 1. Core Partioning
            [DSER, IR_Late] = iceberg_core(testCase.IrObj);
            testCase.verifyEqual(DSER.nChannels, 4);
            testCase.verifyEqual(IR_Late.nChannels, 4);
            testCase.verifyNotEmpty(DSER.time);
            testCase.verifyNotEmpty(IR_Late.time);
            
            % 2. Late rendering Ambisonics
            iAngle = 45; % target azimuthal panning target
            
            % Convert configurations dynamically as done inside `iceberg.m`
            [C,ia,ic] = unique(testCase.ConfigSetup.lsArray);
            ambSetup.ls_dir = [(C(1:end))', zeros(length(C),1)];
            
            signal_amb = iceberg_set_amb(testCase.SignalObj, IR_Late, ambSetup);
            testCase.verifyEqual(signal_amb.nChannels, 24, 'Output from ambisonics should match array size.');
            testCase.verifyNotEmpty(signal_amb.time);
            
            % 3. Early Rendering VBAP
            signal_vbap = iceberg_set_vbap(testCase.SignalObj, DSER, repmat(iAngle, 1, 4), testCase.ConfigSetup);
            testCase.verifyEqual(signal_vbap.nChannels, 4, 'Output from VBAP should match active node setup.');
            
            % 4. Merging Back to physical Matrix array layout (iceberg_merge)
            signal_final = iceberg_merge(testCase.SignalObj, DSER, IR_Late, -40, repmat(iAngle, 1, 4), testCase.ConfigSetup);
            testCase.verifyEqual(signal_final.nChannels, 24, 'Final merged output must match full array configuration.');
            testCase.verifyEqual(length(signal_final.time), length(signal_amb.time), 'Final stream duration should match long Amb array.');
        end
        
    end
end

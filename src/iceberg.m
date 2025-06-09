function iceberg_signal = iceberg(signal, IR, Selected_Angle, configSetup)
% ICEBERG Auralization Method (VBAP + Ambisonics)
% Combines VBAP for direct sound and early reflections with Ambisonics for
% late reverberation in spatial audio rendering.
%
% Inputs:
%   signal          - Input audio signal (itaAudio object)
%   IR              - Impulse responses
%   selectedAngle   - Selected source angle
%   configSetup     - Configuration setup structure
%
% Output:
%   iceberg_signal  - Processed spatial audio signal
[dser, lr] = iceberg_core(IR);

vbapDser         = iceberg_set_vbap(signal, dser, Selected_Angle, configSetup);
% vbapDser         = calibrate_vbap(vbapDser, dser, configurationSetup);

ambisonicsLr   = iceberg_set_amb(signal, lr, configSetup);
% Ambisonics_ERLR = calibrate_ambisonics(Ambisonics_ERLR,level,iAngles,configurationSetup)

vbapAudio = createAudioObject(signal.samplingRate);
ambAudio  = createAudioObject(signal.samplingRate);


activeLsNumbers = configSetup.activeLSNumbers;

% Distribute signals to appropriate loudspeaker channels
for i = 1:length(activeLsNumbers)
    channel = activeLsNumbers(i);
    vbapAudio.time(:, channel) = vbapDser.time(:, i);
    ambAudio.time(:, channel)  = ambisonicsLr.time(:, i);
end

iceberg_signal = ita_add(vbapAudio, ambAudio);

%Pad array with zeros
if length(configSetup.lsArray) > iceberg_signal.dimensions
    iceberg_signal.time(:,iceberg_signal.dimensions+1:length(configSetup.lsArray))...
        = zeros(iceberg_signal.nSamples,...
        length(configSetup.lsArray)-(iceberg_signal.dimensions));
end

%% Helper Functions
    function audioObj = createAudioObject(samplingRate)
        % Creates and initializes an itaAudio object
        audioObj = itaAudio();
        audioObj.samplingRate = samplingRate;
    end
end
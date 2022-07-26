%% It is necessary to set the input output inside ita toolbox
% Please run ita_preferences and adjust in I/O settings the recording and
% the playing device to run trhough portaudio (http://www.playrec.co.uk/)

% Auralization (final 24 channels with content to 4 LS):

Hibrid_4LS

%% Play the resulting file

playPlain(Signal_To_Run)

%% Play function
function playPlain(signal)
[~, playDeviceInfo] = ita_portaudio_deviceID2string(ita_preferences('playDeviceID'));
[~, recDeviceInfo] = ita_portaudio_deviceID2string(ita_preferences('recDeviceID'));
playDeviceID = playDeviceInfo.deviceID;
recDeviceID = recDeviceInfo.deviceID;
iFs = signal.samplingRate;                %Sample Frequency
signal_run = signal;
nCh = 24;
if playrec('isInitialised')==1
    playrec('reset')
    playrec('init', iFs, playDeviceID, recDeviceID)
    playrec('play',signal_run.time,1:nCh);
else
    playrec('init', iFs, playDeviceID, recDeviceID)
    playrec('play',signal_run.time,1:nCh);
end
end
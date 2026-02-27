% angularSpacing = 360/numberOfSpeakers;
% lsAnglesPart1 = mod(iceberglAngles(1):angularSpacing:360, 360);  % 180° to 360° (0°)
% lsAnglesPart2 = angularSpacing:angularSpacing:iceberglAngles(1)-angularSpacing;  % 15° to 165°
% configurationSetup.LSArray = [lsAnglesPart1, lsAnglesPart2]';
% activeLSNumbers = zeros(1,length(configurationSetup.ls_dir));
% for indexx= 1:length(configurationSetup.ls_dir)
%     iArrayP = find(configurationSetup.LSArray==configurationSetup.ls_dir(indexx,1),1);
%     activeLSNumbers(indexx) = iArrayP;
% end
% configurationSetup.activeLSNumbers = activeLSNumbers;
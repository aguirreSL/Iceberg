function setToolboxes()
try
    test = itaAudio;
    clear test
catch
    disp('Indeed we need ITA-Toolbox. We shall install it, dear')
cd ..
if ispc
run([pwd '\Toolboxes\ITA-Toolbox\ita_toolbox_setup.m'])
else
    !git clone https://git.rwth-aachen.de/ita/toolbox.git ../Toolboxes/ITA-Toolbox
    run([pwd '/Toolboxes/ITA-Toolbox/ita_toolbox_setup.m'])
end

cd Calibration
end

try
    ls_dirs = [30 -30 0 110 -110]; % define a 2D 5.0 setup in degrees
    ls_groups = findLsPairs(ls_dirs);
catch e
%     disp(e.message)
   !git clone https://github.com/polarch/Vector-Base-Amplitude-Panning.git ./Toolboxes/polarch/vbap
end

try
    [~, ls_dirs4_rad] = platonicSolid('tetra');
catch e
%     disp(e.message)
    clear ls_dirs
   !git clone https://github.com/polarch/Higher-Order-Ambisonics.git ./Toolboxes/polarch/hoa
end
    clear ls_dirs ls_dirs4_rad ls_groups e
%% Could'nt find VADSHON repository. Downloaded from Imperial College
%add to path
addpath(genpath(pwd))

disp('Required toolboxes installed')

end

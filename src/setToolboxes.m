function setToolboxes()
% SETTOOLBOXES - Checks and installs required audio processing toolboxes
%   This function ensures the ITA-Toolbox, VBAP, and HOA toolboxes are available
%   - Installs missing toolboxes from official repositories
%   - Sets up paths for cloned repositories
%   - Provides clear status messages during execution

%% Check and install ITA-Toolbox
if ~itaToolboxInstalled()
    disp('ITA-Toolbox not found. Installing from RWTH Aachen repository...');
    
    % Clone repository
    itaRepo = 'https://git.rwth-aachen.de/ita/toolbox.git';
    itaPath = fullfile(pwd, 'Toolboxes', 'ITA-Toolbox');
    cloneRepository(itaRepo, itaPath);
    
    % Run setup
    setupFile = fullfile(itaPath, 'ita_toolbox_setup.m');
    if isfile(setupFile)
        run(setupFile);
        disp('ITA-Toolbox installed successfully.');
    else
        error('ITA-Toolbox setup file not found at: %s', setupFile);
    end
end

%% Check and install VBAP toolbox
if ~vbapToolboxInstalled()
    disp('VBAP toolbox not found. Installing from GitHub...');
    vbapRepo = 'https://github.com/polarch/Vector-Base-Amplitude-Panning.git';
    vbapPath = fullfile(pwd, 'Toolboxes', 'polarch', 'vbap');
    cloneRepository(vbapRepo, vbapPath);
    addpath(genpath(vbapPath));
    disp('VBAP toolbox installed successfully.');
end

%% Check and install HOA toolbox
if ~hoaToolboxInstalled()
    disp('HOA toolbox not found. Installing from GitHub...');
    hoaRepo = 'https://github.com/polarch/Higher-Order-Ambisonics.git';
    hoaPath = fullfile(pwd, 'Toolboxes', 'polarch', 'hoa');
    cloneRepository(hoaRepo, hoaPath);
    addpath(genpath(hoaPath));
    disp('HOA toolbox installed successfully.');
end

%% Final path update and cleanup
addpath(genpath(pwd));
disp('All required toolboxes are installed and available.');

% ----------------- Helper Functions -----------------

function installed = itaToolboxInstalled()
    try
        % Minimal check for ITA-Toolbox presence
        itaAudio();
        installed = true;
    catch
        installed = false;
    end
end

function installed = vbapToolboxInstalled()
    % Check for VBAP function
    installed = exist('findLsPairs', 'file') == 2;
end

function installed = hoaToolboxInstalled()
    % Check for HOA function
    installed = exist('platonicSolid', 'file') == 2;
end

function cloneRepository(repoUrl, targetPath)
    % Create parent directories if needed
    [parentPath, ~] = fileparts(targetPath);
    if ~isfolder(parentPath)
        mkdir(parentPath);
    end
    
    % Clone using system-independent command
    if ispc
        system(sprintf('git clone %s "%s"', repoUrl, targetPath));
    else
        system(sprintf('git clone %s ''%s''', repoUrl, targetPath));
    end
    
    % Verify successful clone
    if ~isfolder(targetPath)
        error('Failed to clone repository to: %s', targetPath);
    end
end

end
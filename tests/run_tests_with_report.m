% RUN_TESTS_WITH_REPORT Runs the Iceberg testing suite and generates an HTML report

% Setup testing architecture
import matlab.unittest.TestRunner;
import matlab.unittest.TestSuite;
import matlab.unittest.plugins.TestReportPlugin;
import matlab.unittest.plugins.DiagnosticsRecordingPlugin;

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));

suite1 = TestSuite.fromFile(fullfile(fileparts(mfilename('fullpath')), 'test_native_replacements.m'));
suite2 = TestSuite.fromFile(fullfile(fileparts(mfilename('fullpath')), 'test_iceberg_integration.m'));
suite = [suite1, suite2];

% Create runner
runner = TestRunner.withTextOutput;

% Add plugins for diagnostics and HTML reporting
runner.addPlugin(DiagnosticsRecordingPlugin);

% Create dedicated reports directory
reportFolder = fullfile(fileparts(mfilename('fullpath')), 'reports');
if ~exist(reportFolder, 'dir')
    mkdir(reportFolder);
end

% Set HTML plugin parameters
htmlFile = fullfile(reportFolder, 'Iceberg_Test_Report.html');
plugin = TestReportPlugin.producingHTML(reportFolder, ...
    'MainFile', 'Iceberg_Test_Report.html', ...
    'IncludingCommandWindowText', true, ...
    'IncludingPassingDiagnostics', true);

runner.addPlugin(plugin);

disp('---------------------------------------------------');
disp('Running Iceberg Native Replacements Test Suite...');
disp('---------------------------------------------------');

% Execute tests
result = runner.run(suite);

disp('---------------------------------------------------');
fprintf('HTML Report successfully generated at:\n %s\n', htmlFile);
disp('---------------------------------------------------');

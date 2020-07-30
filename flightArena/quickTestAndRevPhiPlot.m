function quickTestAndRevPhiPlot(relDir)

% this function uses the experiment directory to plot both the
% quickTestTrialPlot and the quickRevPhiTrailPlot
%
% relDir -  Directory contaning the experimental Log file and the Trials Log
% files (if there is more than one trial)

close all

tempDir = dir(relDir);

testInd = find(cellfun(@(x) contains(x, 'trial'), {tempDir.name}));

if isempty(testInd)
    fprintf('No test data in directory\n')
else
    for ii=1:length(testInd)
        
        tempDir2 = dir(fullfile(relDir, tempDir(testInd(ii)).name)); 
        logInd = find(cellfun(@(x) contains(x, 'G4_TDMS_Logs'), {tempDir2.name}));
        
        assert(length(logInd)==1, 'problem with test log file')
        
        logFileName = fullfile(relDir, tempDir(testInd(ii)).name, tempDir2(logInd).name);
        load(logFileName, 'Log')
        
        quickTestTrialPlot(Log)
        
    end
end


eLogInd = cellfun(@(x) contains(x, 'G4_TDMS_Logs'), {tempDir.name});

switch sum(eLogInd)
    case 0 
        error('No experimental file in directory')
    case 1
        eLogFileName = fullfile(relDir, tempDir(eLogInd).name);
        load(eLogFileName, 'Log')
        
        quickRevPhiTrialPlot(Log)
    otherwise
        error('problem with experiment log file')
end


% to flip figures (early trials on top
fh = get(0, 'children');

set(0, 'children', flipud(fh))



end


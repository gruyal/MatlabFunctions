function singleBarSt = generateAlignedSingleBarStDiffWandVwTable(pStruct, pathFlag)

% singleBarSt = generateAlignedSingleBarStwMinMaxDiffWandV(pStruct)
%
% This function take the original protocolStruct from a SingleBarDiagCorr protocol 
% after appear and disappear have been added to the gratingTable and aligns
% it to bar appearance. Resulting structure is organized in position X stimDur structure
% Function also adds baseline subtracted data (baseline calculted by
% prestimBaseWin
%
% This is a simplified version of
% generateAlignedSingleBarStwMinMaxDiffWandV that doesnt generate max and
% min for each position or the norm versions of those but does align
% everything using the table
%
%
% INPUT
%
% pStruct -         protocolStruct for singleBar experiment after relevant
%                   varibles have been added to gratingTable (i.e 'appear' - for frame in whcih bar appears)
% pathFlag -        logical (default 0). marks whether it the recording is from the
%                   OFF (0) or ON (1) pathway. This will determine time
%                   windows for the calculation of max/min responses (since
%                   ON stim in OFF pathway needs longer window
%
% OUTPUT
%
% singleBarSt -         structure with the following fields:
%   .data -             output from alignProtocolDataByTable
%   .subData -          baseline subtracted data. contains following fields:
%       .baseSub -      baseline subtracted data (time and resp)
%       .baseMed -      same as above only median response and not mean
%                       (contains only resp)
%       .baseline -     preStim vector from which baseline was calculted
%       .zeroInd -      index at which point zero time is reached (bar appears)
%       .length -       length of the response vector
%   .resp -             calculted response (could change in future). Based on
%                       sdFac SDs from the baseline (for max in a stimDur
%                       dependent manner). Contains the following fields:

if nargin < 2
    pathFlag = 0;
end

assert(ismember(pathFlag, [0,1]), 'pathFlag should be a logical 0 for T5 and 1 for T4');

preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation
inhibWinFac = 4; % how many time the E window size for detecting I(in old function was to the end)

%since T5 has a more noisy baseline 
if pathFlag
    sdFac = 3;
else
    sdFac = 2.5;
end

sdFacWidMod = 0.2; %factor to be multiflied by floor(width/2) and added to sdFac

postStimWin = [75, 150]; %window in which max/min resp will be calculated (actual window will be 0 to postStimWin+stimdur(in ms) 
maxQ = 0.999;
minQ = 0.001;
smWin = 1000; % smoothing window (to estimate beginning and end of rise/decay phases
respTimeMin = 20; %in ms. if response if found to start before this time it overwrite it 
sampToMsFac = 20; %since data was collcted @ 20KHz

% checking inputs
relTab = pStruct.gratingTable;

assert(ismember('appear', relTab.Properties.VariableNames), ...
           'gratingTable is missing appear variable')

assert(ismember('position', relTab.Properties.VariableNames), 'This function is designed for singleBar protocols only')

alignSt = alignProtocolDataByTable(pStruct, 'appear');

singleBarSt = struct; 

for ii=1:length(alignSt)
    
    singleBarSt.result(ii).data = alignSt(ii); 
    
    tempT = pStruct.gratingTable(ii, :);
    
    assert(isequal(alignSt(ii).table, tempT), 'tables from alignSt and original struct are not equal')
    
    relDat = singleBarSt.result(ii).data.align.mean;
    baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
    zeroInd = find(relDat(:,1) > 0, 1, 'first');

    baseVal = mean(relDat(baseInds(1):baseInds(2), 2));
    baseSubResp = relDat;
    baseSubResp(:,2) = baseSubResp(:,2) - baseVal;

    singleBarSt.result(ii).subData.baseSub = baseSubResp;
    singleBarSt.result(ii).subData.baseline = baseVal;
    singleBarSt.result(ii).subData.zeroInd = zeroInd;
    singleBarSt.result(ii).subData.length = size(baseSubResp, 1);
    
end

singleBarSt.table = pStruct.gratingTable; 



end


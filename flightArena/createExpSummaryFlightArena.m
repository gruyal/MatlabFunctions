function resultSt = createExpSummaryFlightArena(relResDir, plotReportFlag)

% function createExpSummaryFlightArena(relResDir)
%
% This function is designed to collect data from all the flies in the same
% experiment and present a summerized version (TBD)
%
% INPUT
% relResDir  -      relevant results directory which contains all the fly
%                   directories 
% plotReportFlag-   (optional, defualt T) logical. If true plots fixation
%                   and optomotor for each fly (in separate figure)
%
% OUTPUT
% 
% resultSt -        structure contianing the fields
%   .name -         fly name (string)
%   .fixation -     position histograms from intertrial data
%   .time -         timestamp vector
%   .dataCh -       data recorded from the relevant channels (given by
%                   relCh)
%   .length -       length of data vectors
%   .dataChName -   channel names
%   .LmRSD1/2 -     standard deviation for LmR of individual repeats in the
%                   first and second half of a repeat trial (used to
%                   determine whether some repeats should be excluded)
%   .excludedReps - repeats that are excluded
%   .means -        means per stimulus after exclusion 

if nargin < 2
    plotReportFlag = 1; 
end


topDir = dir(relResDir);

relDirInd = find(cellfun(@(x) contains(x, 'fly'), {topDir.name}));

resFileName = '*G4_Processed*'; 
resultSt = struct; 
relCh = 1:5; 
sumCh = [2,3]; % to generate LpR
sdCutoff = 0.1; % determined empirically
optoStim = [21, 23]; %[9,10];  experiment dependent
yyLim = [-6,6];
posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, -0.02, 0.04, 3);
relCol = cbrewer('qual', 'Set1', 3);
arenaXSize = 192;
xPosFix = [1:arenaXSize/2, NaN, (arenaXSize/2+1:arenaXSize) - arenaXSize];


for dd=1:length(relDirInd)
    
    relFN = dir(fullfile(relResDir, topDir(relDirInd(dd)).name, resFileName));
    
    if isempty(relFN)
        warning('directory %s does not contain processed file', topDir(relDirInd(dd)).name)
        continue
    elseif length(relFN) > 1
        warning('directory %s does contains more than one processed file', topDir(relDirInd(dd)).name)
        continue
    else
        
        load(fullfile(relResDir, topDir(relDirInd(dd)).name, relFN.name), 'timeseries', 'timestamps', 'interhistogram', 'histograms', 'channelNames')
        
        datLength = length(timestamps); 
        numStim = size(timeseries,2);
        
        resultSt(dd).name = topDir(relDirInd(dd)).name;
        resultSt(dd).fixation = interhistogram(:, 1:arenaXSize); % since sometimes it creates a longer vector by mistake 
        resultSt(dd).time = timestamps; 
        resultSt(dd).length = datLength;
        resultSt(dd).dataCh = timeseries(relCh, :, :, :); 
        resultSt(dd).LpR = timeseries(sumCh(1), :, :, :) + timeseries(sumCh(2), :, :, :); 
        resultSt(dd).dataChName = channelNames.timeseries(relCh);
        resultSt(dd).LmRSD1 = nanstd(squeeze(timeseries(relCh(1), :, :, 1:floor(datLength/2))), 0, 3); 
        resultSt(dd).LmRSD2 = nanstd(squeeze(timeseries(relCh(1), :, :, floor(datLength/2):end)), 0, 3); 
        resultSt(dd).excludedReps = resultSt(dd).LmRSD1 < sdCutoff | resultSt(dd).LmRSD2 < sdCutoff; % since sometimes flies start/stop half way through
        
        datMeans = nan(length(relCh), numStim, datLength);
        LpRMeans = nan(numStim, datLength);
%         datMedians = datMeans;
        
        for cc=1:length(relCh)
            
            for stim=1:numStim
            
                tempDat = squeeze(timeseries(relCh(cc), stim, ~resultSt(dd).excludedReps(stim, :), :));
                datMeans(cc, stim, :) = nanmean(tempDat); 
%                 datMedians(cc, stim, :) = nanmedian(tempDat); 
                
            end
            
        end
        
        for stim=1:numStim
            
            tempLpR = squeeze(timeseries(sumCh(1), stim, ~resultSt(dd).excludedReps(stim, :), :)) + squeeze(timeseries(sumCh(2), stim, ~resultSt(dd).excludedReps(stim, :), :));
            LpRMeans(stim, :) = nanmean(tempLpR); 
                
        end
        
        resultSt(dd).means = datMeans; 
        resultSt(dd).LpRmeans = LpRMeans; 
%         resultSt(dd).medians = datMedians; 
        
        if plotReportFlag
            figure('name', resultSt(dd).name)
            axh(1) = axes('position', posCell{1});
            hold on 
            fixRes = resultSt(dd).fixation;
            fixRes1 = sum(fixRes(1:floor(size(fixRes,1)/2), :)); % adding a nan in the middle to break the line
            nFixRes = nan(1, arenaXSize+1);
            nFRInd = setdiff(1:arenaXSize+1, arenaXSize/2+1);
            nFixRes(nFRInd) = fixRes1;
            
            plot(xPosFix, nFixRes, 'linewidth', 2, 'color', relCol(1,:))
            
            fixRes2 = sum(fixRes(floor(size(fixRes,1)/2):end, :)); % adding a nan in the middle to break the line
            nFixRes = nan(1, arenaXSize+1);
            nFixRes(nFRInd) = fixRes2;
            plot(xPosFix, nFixRes, 'linewidth', 2, 'color', relCol(2,:))
            hold off
            
            axh(1).XLim = [-arenaXSize/2, arenaXSize/2];
            title('Frame position - fixation only')
            
            
            axh(2) = axes('position', posCell{2});
            tempLMR = squeeze(resultSt(dd).dataCh(relCh(1), setdiff(1:numStim, optoStim), :, :));
            tempLMR = tempLMR(:);
            tempLMR(abs(tempLMR) < 0.1) = nan; 
            histogram(tempLMR)
            title('L minus R histogram - non opto stim')
            
            axh(3) = axes('position', posCell{3});
            hold on 
            line([0,4.5], [0,0], 'linewidth', 1, 'color', 'k')
            plot(timestamps, squeeze(resultSt(dd).means(1,optoStim(1), :)), 'linewidth', 2, 'color', relCol(1,:))
            plot(timestamps, squeeze(resultSt(dd).means(1,optoStim(2), :)), 'linewidth', 2, 'color', relCol(2,:))
            
            hold off
            title('L minus R - optomotor stim')
            
            axh(3).YLim = yyLim;
            
        end
            
        
        
    end
    
end
    
    



end
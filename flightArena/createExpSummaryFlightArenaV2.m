function resultSt = createExpSummaryFlightArenaV2(relResDir, respWin)

% function createExpSummaryFlightArena(relResDir)
%
% This function is a version of createExpSummaryFlightArena since processedFile configuraton changed a bit. 
% This function is designed to collect data from all the flies in the same
% experiment and present a summerized version.
% expanded the output to include normlized lmr and freq
%
% INPUT
% relResDir  -      relevant results directory which contains all the fly
%                   directories 
% respWin -         (optional) 1X2 vector for start and stop of time to analyze
%                   responses over (in seconds) default is [1,2]
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

if nargin<2
    respWin = [1,2]; % in secs
end
    


topDir = dir(relResDir);
preDirInd = cellfun(@(x) ~contains(x, '.'), {topDir.name});
topDir = topDir(preDirInd);

relDirInd = find([topDir.isdir]);

resFileName = 'ProcessedData.mat'; 
metFileName = 'metadata.mat';
resultSt = struct; 
relNames = {'LmR', 'LpR', 'F_chan'}; 
relChInd = zeros(size(relNames));

sdCutoff = 0.1; % determined empirically

relQ = 0.95; 

arenaXSize = 192;


count = 0; 

for dd=1:length(relDirInd)
    
    relFN = dir(fullfile(relResDir, topDir(relDirInd(dd)).name, resFileName));
    metFN = dir(fullfile(relResDir, topDir(relDirInd(dd)).name, metFileName));
    
    
    if isempty(relFN)
        warning('directory %s does not contain processed file', topDir(relDirInd(dd)).name)
        continue
    elseif length(relFN) > 1
        warning('directory %s contains more than one processed file', topDir(relDirInd(dd)).name)
        continue
    else
        
        load(fullfile(relResDir, topDir(relDirInd(dd)).name, relFN.name), 'timeseries', 'timestamps', 'interhistogram', 'channelNames')
        load(fullfile(relResDir, topDir(relDirInd(dd)).name, metFN.name), 'metadata')
        count=count+1;
        
        for nn=1:length(relNames)
            relChInd(nn) = find(strcmp(relNames{nn}, channelNames.timeseries));
        end
        
        datLength = length(timestamps); 
        numStim = size(timeseries,2);
        
        resultSt(count).name = topDir(relDirInd(dd)).name;
        resultSt(count).fixation = interhistogram(:, 1:arenaXSize); % since sometimes it creates a longer vector by mistake 
        resultSt(count).time = timestamps; 
        resultSt(count).lightCyc = metadata.light_cycle;
        resultSt(count).length = datLength;
        resultSt(count).dataCh = timeseries(relChInd, :, :, :); 
        resultSt(count).dataChName = channelNames.timeseries(relChInd);
        resultSt(count).LmRSD1 = nanstd(squeeze(timeseries(relChInd(1), :, :, 1:floor(datLength/2))), 0, 3); 
        resultSt(count).LmRSD2 = nanstd(squeeze(timeseries(relChInd(1), :, :, floor(datLength/2):end)), 0, 3); 
        resultSt(count).excludedReps = resultSt(count).LmRSD1 < sdCutoff | resultSt(count).LmRSD2 < sdCutoff; % since sometimes flies start/stop half way through
        
        datMeans = nan(length(relChInd), numStim, datLength);
        
        for cc=1:length(relChInd)
            
            for stim=1:numStim
            
                tempDat = squeeze(timeseries(relChInd(cc), stim, ~resultSt(count).excludedReps(stim, :), :));
                datMeans(cc, stim, :) = nanmean(tempDat); 
                
            end
            
        end
        
        resultSt(count).means = datMeans; 
        
        
        tempTimeI = arrayfun(@(x) find(resultSt(count).time > x, 1, 'first'), respWin);
        relResp = squeeze(resultSt(count).means(1, :, tempTimeI(1):tempTimeI(2)));
        meanLmR = mean(relResp,2); 
        normFac = quantile(abs(relResp(:)), relQ);
        resultSt(count).meanLmR = meanLmR; 
        resultSt(count).normFacLmR = squeeze(resultSt(count).means(1, :, :)) ./ normFac; 
        normFacResp = resultSt(count).normFacLmR(:, tempTimeI(1):tempTimeI(2));
        resultSt(count).normFacLmRmean = mean(normFacResp, 2); 
        resultSt(count).meanFreq = nanmean(squeeze(resultSt(count).means(3, :, :)),2);
    end
    
end
    
    



end
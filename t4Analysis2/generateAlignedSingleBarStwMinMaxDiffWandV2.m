function [singleBarSt, varargout] = generateAlignedSingleBarStwMinMaxDiffWandV2(pStruct, pathFlag)

% [singleBarSt, varargout] = generateAlignedSingleBarStwMinMaxDiffWandV2(pStruct)
%
% This function take the original protocolStruct from a SingleBarDiagCorr protocol 
% after appear and disappear have been added to the gratingTable and aligns
% it to bar appearance. Resulting structure is organized in position X stimDur structure
% Function also adds baseline subtracted data (baseline calculted by
% prestimBaseWin
%
% Note! this should be a simplified version of generateAlignedSingleBarSt
%
% Note! this function is a modified version of
% generateAlignedSingleBarStwMinMax and was designed to deal also with
% singleBar data that has different bar widths and different bar values
% (also changes sdFac based on width - since expecting stronger resp)
%
% NOTE! changed the time window to detect inhibition since in the
% singleBarDiffW protocol also increased the ISI
%
% NOTE! (2019 02) This function is a modification of generateAlignedSingleBarStwMinMaxDiffWandV
% it includes several differences: 
% 
%   (1) mean was taken of the baseline and a guassian was used to estimate the
%       distribution of the mean - This mainly clalculates the between stim
%       fluctuations. 
%       In this function the SD of each baseline is taken and therefore it
%       estimates within stim fluctuations better
%   (2) the function uses alignProtocolDataByTable2 (which uses
%       getAlignStimDataByTable2). These use gaussian estimation of both
%       pre and post stimulus baseline as a yard stick for excluding noisy
%       repeats
%       getAlignStimDataByTable2 also removes electrical blips
%
%       Also deleted FWHM calculation 
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
%       .max/minVal-    maximum/ minimum value of resp (actually a
%                       quantile max/minQ)
%       .max/minTime -  time at which value is reached (since it is a
%                       quantile and not really max)
%       .maxInd -       index in which max was found
%       .pre/postRespTime - used when fitting exp
%
% baseSt -              (optional) to be used in alignning other protocol after alignning singlebar 
% 
% Note!! 
% additional summary data is added by addNormMaxAndMinToSingleBarAligned to
% the (end+1, end+1) position in the structure


if nargin < 2
    pathFlag = 0;
end

assert(ismember(pathFlag, [0,1]), 'pathFlag should be a logical 0 for T5 and 1 for T4');

preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation
inhibWinFac = 4; % how many time the E window size for detecting I(in old function was to the end)

% Changed the SD calculation - so need to change the threshold 
if pathFlag
    sdFac = 5; 
else
    sdFac = 2.5;
end

sdFacWidMod = 0.2; %factor to be multiflied by floor(width/2) and added to sdFac

% postStimWin = [75, 150]; %window in which max/min resp will be calculated (actual window will be 0 to postStimWin+stimdur(in ms) 
postStimWin = [75, 150] * 2.5; % need to change it for new t4 recording since some cells had slower Taus
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

[alignSt, ~, baseSt] = alignProtocolDataByTable2(pStruct, {'appear', 'appear', 'disappear'});

uPos = unique(relTab.position);
uDur = unique(relTab.stimDur);
uVal = unique(relTab.value);
uWid = unique(relTab.width);

% emptyIndMat = zeros(length(uPos), length(uDur), length(uWid), length(uVal));

singleBarSt = struct;
allBaseline = [];

for ii=1:length(uPos)
    
    for jj=1:length(uDur)
        
        for kk=1:length(uWid)
            
            for mm=1:length(uVal)
                
                relInd = ismember(relTab{:, {'position'; 'stimDur'; 'width'; 'value'}}, [uPos(ii), uDur(jj), uWid(kk), uVal(mm)], 'rows');
                
                if sum(relInd) == 1
                
                    singleBarSt(ii, jj, kk, mm).data = alignSt(relInd);
                    relDat = singleBarSt(ii, jj, kk, mm).data.align.mean;
                    relDatMed = singleBarSt(ii, jj, kk, mm).data.align.median;
                    baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
                    zeroInd = find(relDat(:,1) > 0, 1, 'first');

                    baseVals = relDat(baseInds(1):baseInds(2), 2);
                    baseSubResp = relDat;
                    baseSubResp(:,2) = baseSubResp(:,2) - mean(baseVals); %used to be all baseVals, but this makes more sense
                    baseSubMed = relDatMed(:,2) - mean(baseVals);

                    singleBarSt(ii, jj, kk, mm).subData.baseSub = baseSubResp;
                    singleBarSt(ii, jj, kk, mm).subData.baseSubMed = baseSubMed;
                    singleBarSt(ii, jj, kk, mm).subData.baseline = baseVals;
                    singleBarSt(ii, jj, kk, mm).subData.zeroInd = zeroInd;
                    singleBarSt(ii, jj, kk, mm).subData.length = size(baseSubResp, 1);

                    allBaseline = vertcat(allBaseline, std(baseVals));
                    singleBarSt(ii, jj, kk, mm).empty = 0;
                    
                elseif sum(relInd) == 0 % combination might not exist in presented stimuli
                    
                    singleBarSt(ii, jj, kk, mm).empty = 1;
                    continue
                    
                else
                    error('relInd does not provide a unique combination for a stimulus with pos %d dur %d wid %d and val %d',  ...
                          uPos(ii), uDur(jj), uWid(kk), uVal(mm))
                    
                end
        
            end
            
        end
        
    end
    
end


% checked that allBaseline is not dependent on time (can't do that on
% average - check individual trials and so no connection)

baseSD = quantile(allBaseline, 0.75); % checked 0.5 0.75 and 0.9 and since the distribution has a right tail 0.75 seems a good compromise

% calculting max/min response

minFoundCutoff = 3; 


% run a similar loop to the one below just to find minSDFac that is shared
% for values and durations (can differ for widths)

minSDFacVec = nan(size(uWid)); 
minDurInd =  length(uDur);
minValInd = find(pathFlag == uVal);
minWidInds = find(uWid > 1);

for kk=1:length(minWidInds)
    
    minFound = 0; 
    minSDFac = 0; 

    while minFound < minFoundCutoff
        
        minFound = 0;  % to avoid carrying minFound from previous values
        minSDFac = minSDFac + 0.5; 
        
        if minSDFac > 0.5
            warning('!!!  Increased minSDFac to %3.1f for width %d   !!!', minSDFac, uWid(minWidInds(kk)))
        end

        for ii=1:length(uPos)

            if singleBarSt(ii, minDurInd, minWidInds(kk), minValInd).empty == 1
                continue
            end

            relPostStimW = postStimWin(1+abs(pathFlag - uVal(minValInd))); % if they are the same, index is one o/w 2

            relFac = sdFac + sdFacWidMod * floor(uWid(minWidInds(kk))/2);

            relDat = singleBarSt(ii, minDurInd, minWidInds(kk), minValInd).subData;
            %smoothing the data for a better derivative
            smoothRelDat = smooth(relDat.baseSub(:,2), smWin, 'sgolay', 3);
            smDerDat = diff(smoothRelDat);
            relRespTime = relPostStimW + uDur(minDurInd)*1000; % to convert to ms
            % since this is only for pathFlag == uVal(mm)
            relRespTimeInhib = inhibWinFac*relPostStimW + uDur(minDurInd)*1000; % to convert to ms
                       
            relIndEnd = find(relDat.baseSub(:,1) > relRespTime, 1, 'first');
            relIndEndInhib = find(relDat.baseSub(:,1) > relRespTimeInhib, 1, 'first');
            
            % if data is shorter goes to the end (since these might be isi
            % diff between t5 and t4 new recordings
            if isempty(relIndEndInhib)
                relIndEndInhib = length(relDat.baseSub(:,1));
            end

            relInd = relDat.zeroInd;

            minResp = quantile(relDat.baseSub(relInd:relIndEndInhib,2), minQ);

            if abs(minResp) > baseSD * (relFac - minSDFac) % since inhibition is smaller
                minFound = minFound +1;             
            end

        end

    end
    
    minSDFacVec(minWidInds(kk)) = minSDFac;

end

% if the first position is for width 1, fill the value with the that of the next width 
if isnan(minSDFacVec(1))
    minSDFacVec(1) = minSDFacVec(2); 
end


for mm=1:length(uVal)
    
    for jj=1:length(uDur)
        
        for kk=1:length(uWid)
            
            for ii=1:length(uPos)

                if singleBarSt(ii, jj, kk, mm).empty == 1
                    continue
                end

                relPostStimW = postStimWin(1+abs(pathFlag - uVal(mm))); % if they are the same, index its one o/w 2

                relFac = sdFac + sdFacWidMod * floor(uWid(kk)/2);

                relDat = singleBarSt(ii,jj,kk,mm).subData;
                %smoothing the data for a better derivative
                smoothRelDat = smooth(relDat.baseSub(:,2), smWin, 'sgolay', 3);
                smDerDat = diff(smoothRelDat);
                relRespTime = relPostStimW + uDur(jj)*1000; % to convert to ms
                % change the condition to make it equal for T4 and T5
                if pathFlag == uVal(mm) % pathFlag == 0 && uVal(mm) == 0
                    relRespTimeInhib = inhibWinFac*relPostStimW + uDur(jj)*1000; % to convert to ms
                elseif pathFlag ~= uVal(mm) % pathFlag == 0 && uVal(mm) == 1
                    relRespTimeInhib = relRespTime; % since in T5, with B stim, inhibition only appears when stim is on 
                end
                relIndEnd = find(relDat.baseSub(:,1) > relRespTime, 1, 'first');
                relIndEndInhib = find(relDat.baseSub(:,1) > relRespTimeInhib, 1, 'first');

                relInd = relDat.zeroInd;

                % max caculated by mean response
                maxResp = quantile(relDat.baseSub(relInd:relIndEnd,2), maxQ);
                minResp = quantile(relDat.baseSub(relInd:relIndEndInhib,2), minQ);

                if maxResp > baseSD * relFac
                    singleBarSt(ii,jj,kk,mm).resp.maxVal = maxResp;
                    tempInd = find(relDat.baseSub(relInd:relIndEnd, 2) > maxResp, 1, 'first') + relInd - 1; 
                    singleBarSt(ii,jj,kk,mm).resp.maxTime = relDat.baseSub(tempInd, 1);
                    singleBarSt(ii,jj,kk,mm).resp.maxInd = tempInd;  
                    preMaxInd = find(relDat.baseSub(1:tempInd,2) < maxResp/2, 1, 'last');
                    postMaxInd = find(relDat.baseSub(tempInd:end,2) < maxResp/2, 1, 'first') +tempInd -1;

                    noMax = 0;
                else
                    singleBarSt(ii,jj,kk,mm).resp.maxVal = 0;
                    singleBarSt(ii,jj,kk,mm).resp.maxInd = nan;
                    singleBarSt(ii,jj,kk,mm).resp.maxTime = nan;
                    noMax = 1;
                end

                if singleBarSt(ii,jj,kk,mm).resp.maxTime < respTimeMin
                    singleBarSt(ii,jj,kk,mm).resp.maxVal = 0;
                    singleBarSt(ii,jj,kk,mm).resp.maxInd = nan;
                    singleBarSt(ii,jj,kk,mm).resp.maxTime = NaN;
                    noMax = 1;
                end

                if abs(minResp) > baseSD * (relFac - minSDFacVec(kk)) % since inhibition is smaller
                    singleBarSt(ii,jj,kk,mm).resp.minVal = minResp;
                    tempInd = find(relDat.baseSub(relInd:end, 2) < minResp, 1, 'first') + relInd - 1;
                    singleBarSt(ii,jj,kk,mm).resp.minTime = relDat.baseSub(tempInd, 1);
                    singleBarSt(ii,jj,kk,mm).resp.minInd = tempInd;
                    noMin = 0;
                else
                    singleBarSt(ii,jj,kk,mm).resp.minVal = 0;
                    singleBarSt(ii,jj,kk,mm).resp.minTime = nan;
                    singleBarSt(ii,jj,kk,mm).resp.minInd = nan;
                    noMin = 1;
                end

                if singleBarSt(ii,jj,kk,mm).resp.minTime < 2*respTimeMin % in case min is found before response zero min resp
                    singleBarSt(ii,jj,kk,mm).resp.minVal = 0;
                    singleBarSt(ii,jj,kk,mm).resp.minTime = nan;
                    singleBarSt(ii,jj,kk,mm).resp.minInd = nan;
                    noMin = 1;
                end


                % fill in times for missing min and max for future exp fitting
                if noMax && ~noMin 
                    tempTime = singleBarSt(ii,jj,kk,mm).resp.minTime;
                    tempInd = find(relDat.baseSub(:,1) == tempTime);
                    singleBarSt(ii,jj,kk,mm).resp.preRespTime = NaN; 
                    % estime baginning from value close to baseline
                    altMaxInd = find(relDat.baseSub(1:tempInd, 2) > -1/3*baseSD, 1, 'last');
                    % improve estimated beginning from smooth derivative
                    altMaxInd2 = find(smDerDat(1:altMaxInd) > 0, 1, 'last');
                    singleBarSt(ii,jj,kk,mm).resp.maxTime = relDat.baseSub(altMaxInd2 + floor(smWin/5), 1); %added smWin/5 to correct for smoothing artefact

                    % if the time found for the start moves before the stimulus,
                    % this moves it back to respTimeMin
                    if singleBarSt(ii,jj,kk,mm).resp.maxTime < respTimeMin
                        singleBarSt(ii,jj,kk,mm).resp.maxTime = respTimeMin;
                    end

                    % if time found for beginning is more than the given window for
                    % response, this zeros also the minTime
                    if singleBarSt(ii,jj,kk,mm).resp.maxTime > relRespTime
                        singleBarSt(ii,jj,kk,mm).resp.minVal = 0;
                        singleBarSt(ii,jj,kk,mm).resp.maxVal = 0;
                        singleBarSt(ii,jj,kk,mm).resp.maxInd = nan;
                        singleBarSt(ii,jj,kk,mm).resp.minInd = nan;
                        singleBarSt(ii,jj,kk,mm).resp.minTime = NaN;
                        singleBarSt(ii,jj,kk,mm).resp.maxTime = NaN;
                    end

                elseif ~noMax && noMin 
                    tempTime = singleBarSt(ii,jj,kk,mm).resp.maxTime;
                    tempInd = find(relDat.baseSub(:,1) == tempTime);
                    altMinInd = find(relDat.baseSub(tempInd:end, 2) < 1/3*baseSD, 1, 'first');

                    if isempty(altMinInd)% in case doesn't go down to zero (wasn't replicated for max since is has to start from baseline)
                        altMinInd = length(relDat.baseSub(tempInd:end, 1));
                    end

                    singleBarSt(ii,jj,kk,mm).resp.minTime = relDat.baseSub(altMinInd+tempInd-1, 1);
                    singleBarSt(ii,jj,kk,mm).resp.minInd = altMinInd+tempInd-1;

                elseif noMax && noMin % just for the sake of cleaness

                    singleBarSt(ii,jj,kk,mm).resp.preRespTime = NaN; 
                end

                % filling in relevant pre-response time 
                if ~noMax
                    tempTime = singleBarSt(ii,jj,kk,mm).resp.maxTime;
                    tempInd = find(relDat.baseSub(:,1) == tempTime);
                    preRespInd = find(relDat.baseSub(relInd:tempInd, 2) < 1/3*baseSD, 1, 'last');  % not zero to make sure only the response is captured             

                    if isempty(preRespInd)
                        singleBarSt(ii,jj,kk,mm).resp.preRespTime = respTimeMin;
                    else
                        singleBarSt(ii,jj,kk,mm).resp.preRespTime = relDat.baseSub(preRespInd+relInd-1, 1); 
                    end

                end

            end
            
        end
        
    end
    
end

singleBarSt = addNormMaxAndMinToSingleBarAlignedDiffWandV(singleBarSt, pathFlag);

% adding the values that were used for minSDFac
for ww=1:length(uWid)
    singleBarSt(end, end, ww, end).minSDFac = minSDFacVec(ww); 
end

if nargout == 2
    varargout{1} = baseSt;
end


end


function sigBarSt = addNormMaxAndMinToSingleBarAlignedDiffWandV(alignedSingleBarSt, pathFlag)

% function sigBarSt = addNormMAxAndMinToSingleBarAlignedDiffWandV(alignedSingleBarSt)
%
% This function takes the finished output from generateAlignedSingleBarStwMinMax
% and adds normMax and normMin to it. For all the speed seperately and for
% the total cell. It is used internally in
% generateAlignedSingleBarStwMinMax after the alignement and max
% calculation has been done
%
% Note! this function is a modification of
% addNormMAxAndMinToSingleBarAligned designed to account for singleBar
% stimuli with different widths and values (4D) 
%
% MORE IMPORTANT NOTE! this function does not find the center of E and I
% based on the normalized responses (doesnt even calculate them). Instead
% it calculates a weighted average of maxResp * position 
%
%
% iNPUT
%
% alignedSingleBarSt -      used internally within
%                           generateAlignedSingleBarStwMinMax, but uses the output of that function
%                           as its input
% pathFlag -                inherited from generateAlignedSingleBarStwMinMax
% OUTPUT
%
% sigBarSt -                to the regualr fields that describe each
%                           position, duration combination, this function adds a numPos+1 numDur+1
%                           subStructure that incluedes the following fields:
%   max/min:                extracted max/min from all pos dur combintations 
%   normMax/Min:            normalized max and min across position by
%                           duration
%   maxExt/inhPosInd:       position in which max excitation/inhibition was measured (normalized per speed and summed for designated speeds)
%   maxExt/InhPosVal:       the value of that positon relative to the
%                           center of the stimulus
%   FWHM -                  full width half maximun at each position (if a
%                           response is found)
%
% Note! changed how center is calcultated by moving 0.5/1.5 postion for
% bars width 2/4 (before rounding)



% indices for speeds to indlues in the total normMax/Min computation
% for T5 changed to only the longest duration for both 
% 
% extRespInds = 1:4;
% inhibRespInds = 3:4;

% parameters used for the calculation of the center 
relDur = 0.16; 
relWid = 2;
relWid2 = 4;

if pathFlag
    relVal = 1;  % for T4 
else
    relVal = 0;  % for T5 
end

numObsCutOff = 3; % less then that, uses relWid2 to determine center


sigBarSt = alignedSingleBarSt;

datSiz = size(alignedSingleBarSt);

% I am assuming that only the trialing value dimension will be squeezed and not both value and width 
if length(datSiz) == 3
    datSiz(4) = 1; 
end

allMax = nan(datSiz);
allMin = allMax;
allFWHM = allMax;

allPos = nan(datSiz(1),1);
allDur = nan(datSiz(2),1);
allWid = nan(datSiz(3),1);
allVal= nan(datSiz(4),1);

for ii=1:datSiz(1)
    
    for jj=1:datSiz(2)
        
        for kk=1:datSiz(3)
            
            for mm=1:datSiz(4)
                
                if alignedSingleBarSt(ii,jj,kk,mm).empty == 1
                    continue
                end
        
                allMax(ii,jj,kk,mm) = alignedSingleBarSt(ii,jj,kk,mm).resp.maxVal;
                allMin(ii,jj,kk,mm) = alignedSingleBarSt(ii,jj,kk,mm).resp.minVal;
                if isfield(alignedSingleBarSt(ii,jj,kk,mm).resp, 'FWHM')
                    allFWHM(ii,jj,kk,mm) = alignedSingleBarSt(ii,jj,kk,mm).resp.FWHM;
                end
                
                % so it wont be constantly overwrittten
                if isnan(allPos(ii))
                    allPos(ii) = alignedSingleBarSt(ii,jj,kk,mm).data.table.position; 
                end
                if isnan(allDur(jj))
                    allDur(jj) = alignedSingleBarSt(ii,jj,kk,mm).data.table.stimDur;
                end
                
                if isnan(allWid(kk))
                    allWid(kk) = alignedSingleBarSt(ii,jj,kk,mm).data.table.width;
                end
                
                if isnan(allVal(mm))
                    allVal(mm) = alignedSingleBarSt(ii,jj,kk,mm).data.table.value;
                end
                
            end
            
        end
        
    end
    
end

%use shifted so when taking a weighted mean, position 0 will still be used
% not sure is actually changed anything
shiftedPos = (1:length(allPos))';
shiftBackNum = find(allPos == 0);

relDurInd = find(allDur == relDur);
relWidInd = find(allWid == relWid);
relWidInd2 = find(allWid == relWid2);
relValInd = find(allVal == relVal); 

for kk=1:datSiz(3)
    
    for mm=1:datSiz(4)
        
        tempMax = allMax(:,:,kk,mm);
        tempMin = allMin(:,:,kk,mm);
        
        sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).max = tempMax;
        sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).min = tempMin;
        sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).FWHM = allFWHM(:,:,kk,mm);
%         normMax = tempMax * diag(1./max(tempMax));
%         normMin = tempMin * diag(1./min(tempMin));
%         sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).normMax = normMax;
%         sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).normMin = normMin;
%         maxSum = nansum(normMax(:,end), 2);
%         minSum = nansum(normMin(:,end), 2);
%         [~, maxI] = max(maxSum);
%         [~, minI] = max(minSum);
        
        if allVal(mm) == relVal %find center just based on OFF resp for T5 and ON for T4
            
            wPos = zeros(1, datSiz(2));
            wPosMin = zeros(1, datSiz(2));
            
            for jj=1:datSiz(2)
                
                tempWMax = tempMax(:,jj)'./ nansum(tempMax(:,jj));
                relI = ~isnan(tempWMax);
                wPos(jj) = (tempWMax(relI) * shiftedPos(relI));
                
                tempWMin = tempMin(:,jj)'./ nansum(tempMin(:,jj));
                relIMin = ~isnan(tempWMin);
                wPosMin(jj) = (tempWMin(relIMin) * shiftedPos(relIMin));
                
            end
            
            numObsMax = sum(tempMax > 0);
            numObsMin = sum(tempMin < 0);
            
            sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).wPosMax = wPos - shiftBackNum;
            sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).wPosMin = wPosMin - shiftBackNum;
            sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).wPosMaxNObs = numObsMax;
            sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).wPosMinNObs = numObsMin;
            
        end

        
        
%         sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).maxExtPosInd = maxI;
%         sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).maxExtPosVal = sigBarSt(maxI, 1,kk,mm).data.table.position;
%         sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).minInhPosInd = minI;
%         sigBarSt(datSiz(1)+1, datSiz(2)+1,kk,mm).minInhPosVal = sigBarSt(minI, 1,kk,mm).data.table.position;

    end
    
end

%smarter way would be to look for min observations that are sequential as
%opposed to seperated
if ~isempty(relWidInd2)
    if sigBarSt(end, end, relWidInd, relValInd).wPosMaxNObs(relDurInd) >= numObsCutOff 
%         totMaxPos = round(sigBarSt(end, end, relWidInd, relValInd).wPosMax(relDurInd));
        totMaxPos = round(sigBarSt(end, end, relWidInd, relValInd).wPosMax(relDurInd)-0.5); % -0.5 since the width is 2
    elseif sigBarSt(end, end, relWidInd2, relValInd).wPosMaxNObs(relDurInd) >= numObsCutOff
%         totMaxPos = round(sigBarSt(end, end, relWidInd2, relValInd).wPosMax(relDurInd)) -1;
        totMaxPos = round(sigBarSt(end, end, relWidInd2, relValInd).wPosMax(relDurInd)-1.5);%since it is width 4
    else
        error('cant find wid with enough max observations')
    end
else
%     totMaxPos = round(sigBarSt(end, end, relWidInd, relValInd).wPosMax(relDurInd));
    totMaxPos = round(sigBarSt(end, end, relWidInd, relValInd).wPosMax(relDurInd)-0.5);
end

% encountered this problem in cell 16 (wid2 one pos and one noisy pos
% flipped DS, but wid 4 is very clear)

if ~isempty(relWidInd2) %only if it exists
    if sigBarSt(end, end, relWidInd, relValInd).wPosMinNObs(relDurInd) >= numObsCutOff
%         totMinPos = round(sigBarSt(end, end, relWidInd, relValInd).wPosMin(relDurInd));
        totMinPos = round(sigBarSt(end, end, relWidInd, relValInd).wPosMin(relDurInd)-0.5);
        
        if totMinPos == totMaxPos % if width 2 gives the same minPos and maxPos use width 4 also
            totMinPos = round(sigBarSt(end, end, relWidInd2, relValInd).wPosMin(relDurInd)-1.5);
        end
            
    elseif  sigBarSt(end, end, relWidInd2, relValInd).wPosMinNObs(relDurInd) >= numObsCutOff
%         totMinPos = round(sigBarSt(end, end, relWidInd2, relValInd).wPosMin(relDurInd)) -1;
        totMinPos = round(sigBarSt(end, end, relWidInd2, relValInd).wPosMin(relDurInd)-1.5);
    else
        error('cant find wid with enough min observations')
    end
else
%     totMinPos = round(sigBarSt(end, end, relWidInd, relValInd).wPosMin(relDurInd));
    totMinPos = round(sigBarSt(end, end, relWidInd, relValInd).wPosMin(relDurInd)-0.5);
end

maxPInd = find(allPos == totMaxPos);
minPInd = find(allPos == totMinPos);

sigBarSt(end,end,end,end).maxExtPosInd = maxPInd;
sigBarSt(end,end,end,end).maxExtPosVal = totMaxPos;
sigBarSt(end,end,end,end).minInhPosInd = minPInd;
sigBarSt(end,end,end,end).minInhPosVal = totMinPos; 



end
function newStruct = alignFlightArenaExpStruct(expStruct, verbose)

% This function aligns the OL part of the flight arena experiment using the first pos
% 1 as a common anchor point.
% this function only processes OL data
% 
% added the functionality of checking CL data also, just for marking
% experiments that have no flight in them
%
% INPUT 
%
% expStruct -       Structure, gnerated by parseComplieExpData 
%                   has OL and CL substractures. in OL has .data and
%                   .position and .time (start and stop display) fields (divided by stim X repeat)
% verbose -         logical. if TRUE plots the aligned position functions
%                   for all stim (default 1)
%
% OUTPUT
%
% newStruct -       same dara as original only adding .alignData and
%                   .meanData fields for each simulus. 
%   Note! the above fields are added after the last repeat. E.g if all stim
%   have 3 repeats - align and mean will appear in the 4th positon


freqCh = 5; % used to label trials in which the fly stopped flying
noFlightThresh = 500; % number of samples ( @ 1KHz)

if nargin < 2
    verbose = 1; 
end

posToAlign = 1;
seqLenToFind = 17; % designed to work with the Dec17 exp, where all patterns have 17 frames and all speeds go over all of them

assert(isfield(expStruct, 'dataOL'), 'no Open loop data to align')

olSize = size(expStruct.dataOL);
dSiz = size(expStruct.dataOL(1,1).data,2); % how many channels were recorded

datLenAndZ = nan(olSize(1), olSize(2), 2);

for ii=1:olSize(1)
    
    for jj=1:olSize(2)
        expInd(jj) = ~isempty(expStruct.dataOL(ii,jj).data);
    end
    
    newJJ = find(expInd); % takes out empty expriments
    
    
    for jj=1:length(newJJ)
        
        relData = expStruct.dataOL(ii,jj).data;
        relPos = expStruct.dataOL(ii,jj).position;
        relTime = expStruct.dataOL(ii,jj).time; 
        
        % find the first relevant postion (but a significant apperance of
        % it (since some fixation postions can creep in)
        [spVal, spInd] = SplitVec(relPos(:,2), 'equal', 'firstval', 'first');
        
%         relLens = spLen(spVal == posToAlign);
        relVInds = spInd(spVal == posToAlign);
        
        [lenSpSpVal, indSpSpVal] = SplitVec(spVal, 'consecutive', 'length', 'first');
        tempspvInd = indSpSpVal(find(lenSpSpVal == seqLenToFind, 1, 'first'));
%         [ii,jj]
        assert(spVal(tempspvInd) == 0, 'stimulus beginning was not found - does not start at zero for stim %d repeat %d', ii, jj);
        firstSeqInd = spInd(tempspvInd); % index for the beginning of the first sequecnce whose length is seqLenToFind

        firstInd = find(relVInds > firstSeqInd, 1, 'first');
        alPosInd = relVInds(firstInd);
        
        alPosTime = relPos(alPosInd, 1);
        relTimeInd = find(relData(:,1) > alPosTime, 1, 'first');
        
        expStruct.dataOL(ii,jj).zeroInd = relTimeInd;
        expStruct.dataOL(ii,jj).alignedTime = expStruct.dataOL(ii,jj).time - alPosTime;
        
        datLenAndZ(ii,jj,:) = [size(relData,1), relTimeInd];
        
    end
    
end
        

newLen = round(max(nanmax(datLenAndZ(:, :, 1))) * 1.1); % makes the total length a bit longer to accomadate jitter in timing
lenToCheck = round(newLen/3); % since position is at half sampling and actual stimulus is shorter


for ii=1:olSize(1)
    
    for jj=1:olSize(2)
        expInd(jj) = ~isempty(expStruct.dataOL(ii,jj).data);
    end
    
    newJJ = find(expInd); % takes out empty expriments
    newZero = nanmax(datLenAndZ(ii,:, 2));
    
    alignMat = nan(newLen, length(newJJ), dSiz);
    usefulVec = true(1, length(newJJ));
    
    for jj=1:length(newJJ)
        
        relData = expStruct.dataOL(ii,jj).data;
        relPos = expStruct.dataOL(ii,jj).position;
        
        oldZero = expStruct.dataOL(ii,jj).zeroInd;
        
        for kk=1:size(relData, 2)
            
            respSt.data = relData(:,kk);
            respSt.zeroInd = oldZero;
            
            % checks for stopped flight
            if kk==freqCh
                stopFlightSamp = sum(respSt.data < 1);
                if stopFlightSamp > noFlightThresh
                    usefulVec(jj) = false(1); 
                end
            end
            
            paddedVec = padRespVecGen(respSt, newZero, newLen, 0);
            
            alignData = paddedVec;
            
            % zeros time to the same point
            if kk==1
                alignData = paddedVec - paddedVec(newZero);
                newRelPosT = relPos(:,1) - paddedVec(newZero);
            end
            
            alignMat(:, jj, kk) = alignData; 
            
        end
        
        expStruct.dataOL(ii,olSize(2)+1).alignPos{jj} = [newRelPosT, relPos(:,2)];
    end
    
    % checking aligned positions
    posCheck = zeros(lenToCheck, length(newJJ));
    for jj=1:length(newJJ)
        tempPos = expStruct.dataOL(ii,olSize(2)+1).alignPos{jj};
        tempPosInd = find(tempPos(:,1) > 0, 1, 'first');
        posCheck(:, jj) = tempPos(tempPosInd:tempPosInd+lenToCheck-1, 2);
    end
%     [ii,jj]
    assert(unique(diff(posCheck, 1, 2)) == 0, 'different between position vectors in stimulus %d', ii)    
    
    expStruct.dataOL(ii,olSize(2)+1).usefulFlag = usefulVec;
    expStruct.dataOL(ii,olSize(2)+1).alignZero = newZero;
    expStruct.dataOL(ii,olSize(2)+1).alignData = alignMat;
    
    
    % includes only flight trials in the mean
    useMat = alignMat(:, usefulVec, :);
    expStruct.dataOL(ii,olSize(2)+1).meanData = squeeze(nanmean(useMat, 2));
    expStruct.dataOL(ii,olSize(2)+1).expPos = tempPos(tempPosInd - 5:end, :); % tried 100 at first but keeps going down
    
end
        
       
newStruct = expStruct;

% plots aligned position functions
if verbose
    
    figure
    posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, 0.02, 0.02, [4,3]); % again for the DEc17 experiment 
    axh = gobjects(size(posCell));
    
    for ii=1:olSize(1)
        
        axh(ii) = axes('position', posCell{ii});
        hold on 
        
        tempAlPosCell = expStruct.dataOL(ii, end).alignPos;
        
        cellfun(@(x) plot(x(:,1), x(:,2), 'linewidth', 2), tempAlPosCell)
        
        hold off
        
        axh(ii).YLim = [0, 20];
        
    end
    
end
    
%% CL no flight tagging

clSize = numel(expStruct.dataCL);


for ii=1:clSize
    
    if isempty(expStruct.dataCL(ii).data)
        usefulTag = false(1); 
    else
    
        relDat = expStruct.dataCL(ii).data(:,freqCh);
        stopFlightSamp = sum(relDat < 1);
        if stopFlightSamp > noFlightThresh
            usefulTag = false(1); 
        else
            usefulTag = true(1); 
        end
        
    end
    
    newStruct.dataCL(ii).usefulTag = usefulTag;

end
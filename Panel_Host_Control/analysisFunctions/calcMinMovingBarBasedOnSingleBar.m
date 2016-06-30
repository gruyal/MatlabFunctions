function minMovSt = calcMinMovingBarBasedOnSingleBar(pStruct, sigBarAlignSt)

% function calcMinMovingBarBasedOnSingleBar(pStruct, sigBarAlignSt)
%
% This function uses the singlebar flashes data to reconstruct the
% minMoving bar results based on simple linear summation
%
% INPUT
% pStruct -         protocolStruct for minMovingBar (after appear,disappear
%                   and framePerStep have been added to gratingTable
% sigBarAlignSt -   strcture generated using generateAlignedSingleBarSt
%
% OUTPUT
% TBD

preStimBaseWin = [-250, -50]; % ms before stim to use for baseline calculation

relTab = pStruct.gratingTable;

assert(all(ismember({'appear'; 'disappear'; 'framePerStep'}, relTab.Properties.VariableNames)), ...
           'gratingTable is missing appear variable')

alignSt = alignProtocolDataByTable(pStruct, 'appear');

allPosSig = zeros(1,size(sigBarAlignSt,1));

for ii=1:size(sigBarAlignSt,1)
    allPosSig(ii) = sigBarAlignSt(ii,1).data.table.position;
end

allDurSig = zeros(size(sigBarAlignSt,2), 1);
for ii=1:size(sigBarAlignSt,2)
    allDurSig(ii) = sigBarAlignSt(1, ii).data.table.stimDur;
end


uPair = unique(relTab.pairInd);
uDur = unique(relTab.stepDur);
uDir = unique(relTab.direction);

assert(all(uDur == allDurSig), 'durations are not equal beteween singleBar and minMoving')

minMovSt = struct;

for ii=1:length(uPair)
    
    for jj=1:length(uDur)
        
        for kk=1:length(uDir)
        
            relInd = ismember(relTab{:, {'pairInd'; 'stepDur'; 'direction'}}, [uPair(ii), uDur(jj), uDir(kk)], 'rows');
        
            minMovSt(ii, jj, kk).data = alignSt(relInd);
                
            relDat = minMovSt(ii, jj, kk).data.align.mean;
            baseInds = arrayfun(@(x) find(relDat(:,1) > x, 1, 'first'), preStimBaseWin);
            zeroInd = find(relDat(:,1) > 0, 1, 'first');
                
            baseVal = mean(relDat(baseInds(1):baseInds(2), 2));
            baseSubResp = relDat;
            baseSubResp(:,2) = baseSubResp(:,2) - baseVal;
                
            minMovSt(ii, jj, kk).subData.baseSub = baseSubResp;
            minMovSt(ii, jj, kk).subData.baseline = baseVal;
            minMovSt(ii, jj, kk).subData.zeroInd = zeroInd;
            minMovSt(ii, jj, kk).subData.length = size(baseSubResp, 1);
            
        end
        
    end
    
end

% calculating linear sum

for ii=1:length(uPair)
    
    for jj=1:length(uDur)
        
        startPos = minMovSt(ii,jj,1).data.table.startPos;
        stopPos = minMovSt(ii,jj,1).data.table.stopPos;
        
        relPos = min(startPos, stopPos):max(startPos, stopPos);
        
        relPosInd = ismember(allPosSig, relPos);
        
        if sum(relPosInd) < length(relPos)
            warning('sigBarPos is missing %d positions for pair %d', length(relPos)-sum(relPosInd), ii)
            continue
        end
        
        relSigBarSt = sigBarAlignSt(relPosInd, jj); %since duration was checked to be the same
        
        for kk=1:length(uDir)
            
            relMovSt = minMovSt(ii,jj,kk);
            relMovTable = relMovSt.data.table;
            relMovPos = relMovTable.appear:relMovTable.framePerStep:relMovTable.disappear-1;
            relMovInd = relMovSt.data.align.meanPos(ismember(relMovSt.data.align.meanPos(:,2), relMovPos), 1);
            
            if uDir(kk) == -1
                relMovInd = flipud(relMovInd); 
            end
            
            assert(length(relMovInd) == length(relSigBarSt), 'singlebar positions and minMov positions do not overlap')
            
            totLen = relMovSt.subData.length;
            shiftedVecs = zeros(totLen, length(relMovInd));
            
            for pp=1:length(relMovInd)
                
                shiftedVecs(:,pp) = padRespVec(relSigBarSt(pp), relMovInd(pp), totLen);
                
            end
                
            minMovSt(ii,jj,kk).linSum = [relMovSt.subData.baseSub(:,1), sum(shiftedVecs, 2)];
            
        end
        
    end
    
end




end
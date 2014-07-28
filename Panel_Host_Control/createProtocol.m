function protocolStruct = createProtocol(protocolStruct)

% function protocolStruct = createProtocol(protocolStruct)
%
% Thsi function uses the data within protocolStruct to add a field with the
% data required for presenting the stim to the same structure.
%
% INPUT
%
% protocolStruct -  should include the following fields:
%
% .gratingStruct -  Structure with all the required fields (see
%                   generateGratingFrame) to generate N gratings.
% .masksStruct -     Structure with 2 fields for each mask .type and .radius to generate M masks 
%                   if M is equal to N and interleave is FALSE then each
%                   gratingSeq is presented in each mask according to their order. If M and N
%                   are different and interleave is TRUE all gratingSeq are presented in all
%                   the masks. 
%                   NOTE! if 2 of the 3 (stimSeq, masks, orientations) have only one value in them interleave can still be negative                    
% .orientations -    1XP vector of values 0:1:7 (multiples of 45 degrees that
%                   will be applied to the grating+mask combination. If P is equal or unequal
%                   to N, will be resolved the same way the masks are resolved.
% .maskPositions -   LX2 matrix of positions onto which the center of the
%                   gratingSeq+mask+orientation will be mapped. Position should be given in
%                   arena coordinates (X pixel, Y pixel). All gratingSeq+mask+orientation
%                   combinations will be presented in all the positions given
%                   positions in the order given or randomize them.
% .repeats -         number of times to repeat the whole protocol.
% .interleave -     Flag that indicates how to combine gratingSeqs,
%                   masks and orientations (see above). 
%                   0 -     all should be of equal length and are not
%                           interleaved
%                   1 -     gratingSeq and masks are of equal length, and
%                           combined one to one but orientation is interleaved.
%                   2 -     all is interleaved
%                   3 -     for cases in which maskPositions is a cell
%                           array of trajectories. gratingStruct can contain several gratings, 
%                           but they all should have only one frame. The resulting interleave of 
%                           gratingXmaskXorientation is then moved using maskPosition trajectories
%                   4 -     maskPositions is again a cell array of
%                           trajectories, but gratingSeq, masks and/or orientations
%                           can have the same number of values, and they will be
%                           combined in a one to one correspondance. 
% .randomize -       structure with 4 logical fields 
%                   .gratingSeq 
%                   .masks 
%                   .orientations
%                   .maskPositions. 
%                   Each field is a logical flag for whether or not to
%                   randomize the componentit is describing. If interleave
%                   is TRUE, randomize will be used only after the
%                   approprite mask was applied to a gratingSeq and the
%                   orientation was set
%
% The .randomize structure can have 2 additional fields:
% .randomize.seed - RNG Seed that was used in generating the structure
% .randomize.useSeed - logical. for whether the current function should
%
%                   reuse the seed to generate the current structure
% .intFrames -      number of empty frames (uniform GS level as background)
%                   between stimuli (each combination of grating, mask,
%                   orientation and position)
% .funcHand -       Function handle that is used in interperting the
%                   grating structure (fed into generateGratingBaseSeq)
% .baseParameters    (OPTIONAL). if exist should have .gsLevel and
%                   .background fields (otherwise assume gs3 and background 0.49)
%
% OUTPUT
% for each repeat adds the following fields to the structure 
% .mats -           a cell array of 32X96XL matrices of the stimuli to feed into dumpframe 
%                   (L depends on how gratingSeq, masks, orientations and
%                   mask positions are interleaved)
% .relInds -        LX4 matrix of indices, referencing the relationship between 
%                   the final stimulus matrix and the original order of the inputs. 
%                   1st column is stimSeq index; 2nd is masks; 3rd is
%                   orientation index; and 4th is mask position index
% .randomize.seed - seed value that was used for the random number
%                   generator
% .freqCorr -       1XL vector of frequency correction numbers (if protocol
%                   has several stimuli (gratingXmask combinations) that have different
%                   spatial frequencies, this using correction will allow to present them at the same spatial frequency 
%
% NOTE! In the non-interleave case, the grating is
% combined with a mask and orientation, then mask positions are applied on
% all the combination and then it is randomized if needed. 
% In the interleave case first grating, masks and orientations are combined. 
% They are randomized if needed, mask positions are applied, and then a
% second randomization occurs (if randomize.maskPositions is T)



%checking inputs 

%arenaSiz = [32,96];


protocolStruct = checkProtocolStruct(protocolStruct);

%Generating the grating sequences
stimSeqCell = cell(1, length(protocolStruct.gratingStruct));
tempFreqCorr = zeros(1, length(stimSeqCell));

%finding the width field (since 4Bar uses different width fields)
fNames = fieldnames(protocolStruct.gratingStruct);
fnInd = cellfun(@(x) ~isempty(x), strfind(fNames, 'width'));
widthCell = fNames(fnInd);
if isempty(widthCell)
    warning('No width information in gratingStructure - freqCorr is irrelevant')
end


for ii=1:length(protocolStruct.gratingStruct)
    stimSeqCell{ii} = generateGratingBaseSeq(protocolStruct.gratingStruct(ii), protocolStruct.funcHand);
    tempFreqCorr(ii) = sum(cellfun(@(x) getfield(protocolStruct.gratingStruct(ii), x), widthCell));
end
protocolStruct.stimSeqCell = stimSeqCell;

% This can be used to correct for temporal freqeuncy changes associated
% with spatial freqency changes (e.g. slow down narrower grating frames by
% this amount
tempFreqCorr = max(tempFreqCorr)./tempFreqCorr;


numSeqs = length(protocolStruct.gratingStruct);
numMasks = length(protocolStruct.masksStruct);
numOrt = length(protocolStruct.orientations);
if iscell(protocolStruct.maskPositions)
    numMaskPos = length(protocolStruct.maskPositions);
else
    numMaskPos = size(protocolStruct.maskPositions, 1);
end
numReps = protocolStruct.repeats;

mats = cell(1,numReps);
relInds = mats;
freqCorr = mats;


frameSize = size(stimSeqCell{1}(:,:,1));
masksMat = zeros([frameSize, numMasks]);
% Generating masks
for ii=1:numMasks
    masksMat(:,:,ii) = generateBaseMask(protocolStruct.masksStruct(ii).type, protocolStruct.masksStruct(ii).radius);
end

protocolStruct.masksMat = masksMat;


% This function differentiates between zeros in pattern and zeros in mask
% by setting pattern zeros to be a bit over zero (implemented in
% generateGratingFrame)
if min(cellfun(@(x) min(x(:)), protocolStruct.stimSeqCell)) == 0
    error('stimSeqs minmal value should be above 0 so it will not be set to background levels after rotation')
end

% maybe should move this into checkProtocol
if isfield(protocolStruct, 'baseParameters')
    gsL = protocolStruct.baseParameters.gsLevel;
    bkgdL = protocolStruct.baseParameters.background;
else
    gsL = 3;
    bkgdL = 0.49; %so it would be rounded to 3 in gs3
end

bkgdVal = round(bkgdL * (2^gsL-1));

% Seeding the random number generator so that it can be repeated if needed
% Or using the exsiting seed

if isfield(protocolStruct, 'randomize')
    if isfield(protocolStruct.randomize, 'useSeed')
        if protocolStruct.randomize.useSeed
            rng(protocolStruct.randomize.seed)
        else
            seed = floor(now);
            protocolStruct.randomize.seed = seed;
        end
    else
        seed = floor(now);
        protocolStruct.randomize.seed = seed;
    end
    seed = floor(now);
    protocolStruct.randomize.seed = seed;
end
    



% If gratingSeq, masks, and orientation are to be matched 1 to 1
% or when 2 of the 3 have just one value
switch protocolStruct.interleave
    
    case 0  % DOES NOT INTERLEAVE (ONE TO ONE BETWEEN GRATINGSEQ, MASKS, AND ORIENTATIONS)
       
      
        [maxInd, relInd] = max([numSeqs,numMasks, numOrt]);
        allRand = [protocolStruct.randomize.gratingSeq, protocolStruct.randomize.masks, protocolStruct.randomize.orientations];
        relRand = allRand(relInd);
     
        if numSeqs == 1
             ssI=ones(1, maxInd);
        else
             ssI=1:numSeqs;
        end
     
        if numMasks == 1
         mkI=ones(1, maxInd);
     else
         mkI=1:numMasks;
        end
     
        if numOrt == 1
             orI=ones(1, maxInd);
        else
             orI=1:numOrt;
        end
        
        numOnes = sum([numSeqs, numMasks, numOrt] == 1);
        assert(numOnes == 2 || sum(diff([numSeqs, numMasks, numOrt])) == 0, ...
                'For interleave 0, either gratingSeqs, masks and orientations are of same length, or 2 of 3 are of length 1')
     
     
        % combining gratingSeq with masks and orientations
        rotSeqs = cell(1, maxInd);
        for ii=1:maxInd
             tempSeq = protocolStruct.stimSeqCell{ssI(ii)};
            for jj=1:size(tempSeq,3)
                tempStim = tempSeq(:,:,jj).*protocolStruct.masksMat(:,:,mkI(ii));
                tempStim = imrotate(tempStim, 45*protocolStruct.orientations(orI(ii)), 'nearest', 'crop');
                %sets background to desired level
                tempStim(tempStim == 0) = bkgdVal;
                rotSeqs{ii}(:,:,jj) = round(tempStim);
            end
        end
     
        clear tempStim %otherwise cant be used as cell later on
     
        % mapping to arena and randomizing
        [relRangeX, relRangeY] = getArenaMaskTransform(protocolStruct.maskPositions);
        emptyIntFrames = ones(size(relRangeX, 2), size(relRangeY, 2), protocolStruct.intFrames) * bkgdVal;
        % Crop the full image so that it will appear in desired position
        
        for jj=1:maxInd
            for kk=1:numMaskPos
                relCrds = {relRangeX(kk,:); relRangeY(kk,:)};
                addEmpty = cat(3, emptyIntFrames(:,:,1), rotSeqs{jj}(relCrds{1}, relCrds{2}, :)); %adds one empty frame in the beginning 
                tempStim{jj}{kk} = cat(3, addEmpty, emptyIntFrames); % adds intFrames empty in the end
            end
        end
     
     
        for ii=1:numReps
        
            %randomize if needed
            stimInds = randomizeMatrixInds([maxInd, numMaskPos], ...
                       [relRand, protocolStruct.randomize.maskPositions]);
       
            tempStimCell = arrayfun(@(x) tempStim{stimInds(x,1)}{stimInds(x,2)}, 1:size(stimInds,1), 'uniformoutput', 0);       
            mats{ii} = cellfun(@(x) addBlinker(x, gsL), tempStimCell, 'uniformoutput', 0);
            relInds{ii} = [ssI(stimInds(:,1))', mkI(stimInds(:,1))', orI(stimInds(:,1))', stimInds(:,2)];
            freqCorr{ii} = reshape(tempFreqCorr(ssI(stimInds(:,1))), [],1);
            %reshape makes sure lengths and freqCorr are same shape even
            %when there are multiple inputs
            
        end
        
    case 1  % GRATINGSEQ AND MASKS ARE NOT INTERLEAVED (ONE TO ONE) BUT THE COMBINATIONS ARE INTERLEAVED WITH ALL ORIENTATIONS 
        
        
        
        
        assert(numSeqs == numMasks, 'For interleave 1, gratingSeqs and masks should be of same length')
        
        rotSeqs = cell(numSeqs, numOrt);
        
        for ii=1:numSeqs
            tempSeq = protocolStruct.stimSeqCell{ii};
            for jj=1:size(tempSeq, 3)
                tempStim = tempSeq(:,:,jj).*protocolStruct.masksMat(:,:,ii);
                for kk=1:numOrt
                    tempStim2 = imrotate(tempStim, 45*protocolStruct.orientations(kk), 'nearest', 'crop');
                    %sets background to desired level
                    tempStim2(tempStim2 == 0) = bkgdVal;
                    rotSeqs{ii, kk}(:,:,jj) = round(tempStim2);
                end
            end
        end
        
        
        % mapping to arena and randomizing
        [relRangeX, relRangeY] = getArenaMaskTransform(protocolStruct.maskPositions);
        emptyIntFrames = ones(size(relRangeX, 2), size(relRangeY, 2), protocolStruct.intFrames) * bkgdVal;
        relRand = [protocolStruct.randomize.gratingSeq, protocolStruct.randomize.orientations]; % since gratingSeq and masks are not interleaved
        
        for ii=1:numReps
            
            stimInds = randomizeMatrixInds([numSeqs, numOrt], relRand);
        
            tempStimCell = arrayfun(@(x) rotSeqs{stimInds(x,1), stimInds(x,2)}, 1:size(stimInds,1), 'uniformoutput', 0);
            secStimInd = randomizeMatrixInds([length(tempStimCell), numMaskPos], [1, protocolStruct.randomize.maskPositions]);
            tempStimMat = cell(1, size(secStimInd, 1)); %size of secStimInd is the actual number of individual stimuli after they have been combined
            
            for jj=1:size(secStimInd, 1)
                relCrds = {relRangeX(secStimInd(jj,2),:); relRangeY(secStimInd(jj,2),:)};
                addEmpty = cat(3, emptyIntFrames(:,:,1), tempStimCell{secStimInd(jj,1)}(relCrds{1}, relCrds{2},:));
                stimWBuffer = cat(3, addEmpty, emptyIntFrames);
                tempStimMat{jj} = stimWBuffer;
            
            end
        
            presentedInd = [stimInds(secStimInd(:,1), [1,1,2]), secStimInd(:,2)]; % since gratingSeqs and masks have the same ind
            
            mats{ii} = cellfun(@(x) addBlinker(x, gsL), tempStimMat, 'uniformoutput', 0);
            relInds{ii} = presentedInd;
            freqCorr{ii} = reshape(tempFreqCorr(presentedInd(:,1)), [],1);
            %reshape makes sure lengths and freqCorr are same shape even
            %when there are multiple inputs
        end
        
    case 2  % INTERLEAVE ALL - COMBINING GRATINGSEQ WITH MASKS AND ORIENTATIONS IN ALL COMBINATIONS 
        
        
        rotSeqs = cell(numSeqs, numMasks, numOrt);
        for ii=1:numSeqs
             tempSeq = protocolStruct.stimSeqCell{ii};
            for jj=1:numMasks
                 tempMask = protocolStruct.masksMat(:,:,jj);
                for kk=1:numOrt
                     tempOri = 45*protocolStruct.orientations(kk);
                    for mm=1:size(tempSeq,3)
                        tempStim = tempSeq(:,:,mm).*tempMask;
                        tempStim = imrotate(tempStim, tempOri, 'nearest', 'crop');
                        %  sets background to desired level
                        tempStim(tempStim == 0) = bkgdVal;
                        rotSeqs{ii, jj, kk}(:,:,mm) = round(tempStim);
                    end
                end
            end
        end
    
        %   randomize if needed 
        relRandInd = [protocolStruct.randomize.gratingSeq, protocolStruct.randomize.masks, protocolStruct.randomize.orientations];
        % mapping to arena and randomizing
        [relRangeX, relRangeY] = getArenaMaskTransform(protocolStruct.maskPositions);
        emptyIntFrames = ones(size(relRangeX, 2), size(relRangeY, 2), protocolStruct.intFrames) * bkgdVal;
    
        for ii=1:numReps
            
            stimInds = randomizeMatrixInds2([numSeqs, numMasks, numOrt], relRandInd);
        
            tempStimCell = arrayfun(@(x) rotSeqs{stimInds(x,1), stimInds(x,2), stimInds(x,3)}, 1:size(stimInds,1), 'uniformoutput', 0);
            secStimInd = randomizeMatrixInds([length(tempStimCell), numMaskPos], [1, protocolStruct.randomize.maskPositions]);
            tempStimMat = cell(1, size(secStimInd, 1));  %size of secStimInd is the actual number of individual stimuli after they have been combined
        
            for jj=1:size(secStimInd, 1)
                relCrds = {relRangeX(secStimInd(jj,2),:); relRangeY(secStimInd(jj,2),:)};
                addEmpty = cat(3, emptyIntFrames(:,:,1), tempStimCell{secStimInd(jj,1)}(relCrds{1}, relCrds{2},:));
                stimWBuffer = cat(3, addEmpty, emptyIntFrames);
                tempStimMat{jj} = stimWBuffer;
            
            end
        
            presentedInd = [stimInds(secStimInd(:,1), :), secStimInd(:,2)];
        
            mats{ii} = cellfun(@(x) addBlinker(x, gsL), tempStimMat, 'uniformoutput', 0);
            relInds{ii} = presentedInd;
            freqCorr{ii} = reshape(tempFreqCorr(presentedInd(:,1)), [],1);
            %reshape makes sure lengths and freqCorr are same shape even
            %when there are multiple inputs
        
        end
     
    case 3  % GRATINGSEQ OF LENGTH 1 INTERLEAVED WITH MASKS AND ORIENTAITONS APPLIED TO EVERY MASKPOSITION CELL (TRAJECTORY) 
        
        % make sure each stimSeq contains only one frame
        stimLengths = cellfun(@(x) size(x, 3), protocolStruct.stimSeqCell);
        assert(prod(stimLengths) == 1, 'If interleave equals 3, each gratingSeq should be just one frame')
        
        % make sure maskPositions is a cell array of trajectories
        assert(iscell(protocolStruct.maskPositions), 'If interleave is equal to 3, maskPositions should be a cell array') 
        maskPosLen = cellfun(@(x) size(x, 2), protocolStruct.maskPositions);
        assert(max(maskPosLen) == 2 && min(maskPosLen) == 2, 'maskPositions should be a cell array of NiX2 position')
        
        
        [maxVal, relInd] = max([numSeqs,numMasks, numOrt]);
        allRand = [protocolStruct.randomize.gratingSeq, protocolStruct.randomize.masks, protocolStruct.randomize.orientations];
        relRand = allRand(relInd);
     
        if numSeqs == 1
             ssI=ones(maxVal,1);
        else
             ssI=(1:numSeqs)';
        end
     
        if numMasks == 1
            mkI=ones(maxVal,1);
        else
            mkI=(1:numMasks)';
        end
     
        if numOrt == 1
             orI=ones(maxVal,1);
        else
             orI=(1:numOrt)';
        end
        
        numOnes = sum([numSeqs, numMasks, numOrt] == 1);
        numMax = sum([numSeqs, numMasks, numOrt] == maxVal);
        
        % added this to deal with the case where max is 1 (and therefore
        % condition in assert would be wrong)
        if maxVal == 1
            numMax = 0;
        end
            
        assert(numOnes+numMax == 3, ...
                'For interleave 3, either gratingSeqs, masks and orientations are of same length, or of length 1')
        
        rotSeqs = cell(maxVal, 1);
        for ii=1:maxVal
             tempSeq = protocolStruct.stimSeqCell{ssI(ii)};
             tempMask = protocolStruct.masksMat(:,:,mkI(ii));
             tempOri = 45*protocolStruct.orientations(orI(ii)); 
             tempStim = tempSeq(:,:,1).*tempMask; % since it was already verified that tempSeq is 1 frame long
             tempStim = imrotate(tempStim, tempOri, 'nearest', 'crop');
             %  sets background to desired level
             tempStim(tempStim == 0) = bkgdVal;
             rotSeqs{ii}(:,:,1) = round(tempStim);
        end
    
        % mapping to arena and randomizing
        [relRangeX, relRangeY] = getArenaMaskTransform(protocolStruct.maskPositions);
        emptyIntFrames = ones(size(relRangeX{1}, 2), size(relRangeY{1}, 2), protocolStruct.intFrames) * bkgdVal;
        
        for ii=1:numReps
            
            stimInds = randomizeMatrixInds([maxVal, numMaskPos], [relRand, protocolStruct.randomize.maskPositions]);
            tempStimMat = cell(1, size(stimInds, 1)); %size of stimInd is the actual number of individual stimuli after they have been combined
            
            for jj=1:size(stimInds, 1)
                concatStim = zeros(size(relRangeX{1}, 2), size(relRangeY{1}, 2), ...
                                   size(relRangeX{stimInds(jj, 2)}, 1));
                baseStim = rotSeqs{stimInds(jj, 1)}; % a single frame
                for kk=1:size(relRangeX{stimInds(jj, 2)}, 1)
                    relCrds = {relRangeX{stimInds(jj,2)}(kk,:); relRangeY{stimInds(jj,2)}(kk,:)};
                    concatStim(:,:,kk) = baseStim(relCrds{1}, relCrds{2},1); %should contain only stims of length 1
                end
                addEmpty = cat(3, emptyIntFrames(:,:,1), concatStim);
                stimWBuffer = cat(3, addEmpty, emptyIntFrames);
                tempStimMat{jj} = stimWBuffer;
            end
        
            presentedInd = [ssI(stimInds(:,1)), mkI(stimInds(:,1)), orI(stimInds(:,1)), stimInds(:,2)];
        
            mats{ii} = cellfun(@(x) addBlinker(x, gsL), tempStimMat, 'uniformoutput', 0);
            relInds{ii} = presentedInd;
            freqCorr{ii} = reshape(tempFreqCorr(presentedInd(:,1)), [],1);
            %reshape makes sure lengths and freqCorr are same shape even
            %when there are multiple inputs
        
        end
    
    case 4  % GRATINGSEQ FRAMES, MASKS AND/OR ORIENTATIONS HAVE THE SAME NUMBER OF ELEMENTS AS MASKPOSITION AND ARE APPLIED ONE TO ONE
        
        if iscell(protocolStruct.maskPositions)
            trajecLengths = cellfun(@(x) size(x, 1), protocolStruct.maskPositions);
            gratingLengths = cellfun(@(x) size(x, 3), stimSeqCell);
        else
            error('If interleave Equals 4 maskPositions should be a cell array of trajectories')
        end
        
        assert(length(unique(trajecLengths))  == 1, 'In interleave 4 all trajectories should have the same length')
        assert(length(unique(gratingLengths)) == 1, 'In interleave 4 all gratingSeqs should have the same length')
        
        sharedLen = trajecLengths(1);
        
        if sharedLen == 1
            error('If each maskPositions cells contain one element, do not use interleave 4')
        else
            maskPosI = 1:sharedLen;
        end
        
        if gratingLengths(1) == 1
            gtframeI = ones(1, sharedLen);
        else
            gtframeI = 1:gratingLengths(1);
        end
        
        if numMasks == 1
            mkI=ones(1, sharedLen);
        else
            mkI=1:numMasks;
        end
     
        if numOrt == 1
            orI=ones(1, sharedLen);
        else
            orI=1:numOrt;
        end
        
        ssI = 1:numSeqs;
        
        allLengths = cellfun(@length, {maskPosI, gtframeI, mkI, orI});
        assert(unique(allLengths) == 1, 'In interleave 4, gratingSeq, masks, and orientations should have either one element or the same number as in a maskPosition cell')
        
        % mapping to arena and randomizing
        [relRangeX, relRangeY] = getArenaMaskTransform(protocolStruct.maskPositions);
        emptyIntFrames = ones(size(relRangeX{1}, 2), size(relRangeY{1}, 2), protocolStruct.intFrames) * bkgdVal;
        % Crop the full image so that it will appear in desired position
        
        finSeq = cell(numSeqs, numMaskPos);
        
        for ii=1:numSeqs
            tempSeq = protocolStruct.stimSeqCell{ii};
            for jj=1:numMaskPos
                concatStim = zeros(size(relRangeX{1}, 2), size(relRangeY{1}, 2), sharedLen);
                for kk=1:sharedLen
                    tempStim = tempSeq(:,:,kk).*protocolStruct.maskMat(:,:,mkI(kk));
                    tempStim = imrotate(tempStim, 45*protocolStruct.orientations(orI(kk)), 'nearest', 'crop');
                    %sets background to desired level
                    tempStim(tempStim == 0) = bkgdVal;
                    relCrds = {relRangeX{jj}(kk,:); relRangeY{jj}(kk,:)};
                    concatStim(:,:,kk) = tempStim(relCrds{1}, relCrds{2});
                end    
                addEmpty = cat(3, emptyIntFrames(:,:,1), concatStim);
                stimWBuffer = cat(3, addEmpty, emptyIntFrames);
                finSeq{ii, jj} = stimWBuffer;
            end
        end
        
        [maxInd, relInd] = max([numSeqs,numMasks, numOrt]);
        allRand = [protocolStruct.randomize.gratingSeq, protocolStruct.randomize.masks, protocolStruct.randomize.orientations];
        relRand = allRand(relInd);
        
        % mapping to arena and randomizing
        for ii=1:numReps
        
            %randomize if needed
            stimInds = randomizeMatrixInds([maxInd, numMaskPos], ...
                       [relRand, protocolStruct.randomize.maskPositions]);
       
            tempStimCell = arrayfun(@(x) finSeq{stimInds(x,1),stimInds(x,2)}, 1:size(stimInds,1), 'uniformoutput', 0);       
            mats{ii} = cellfun(@(x) addBlinker(x, gsL), tempStimCell, 'uniformoutput', 0);
            relInds{ii} = [ssI(stimInds(:,1))', mkI(stimInds(:,1))', orI(stimInds(:,1))', stimInds(:,2)];
            freqCorr{ii} = reshape(tempFreqCorr(ssI(stimInds(:,1))), [],1);
            %reshape makes sure lengths and freqCorr are same shape even
            %when there are multiple inputs
            
        end 
        
end

totMats = [mats{:}]';
totRelInds = vertcat(relInds{:});
totFreqCorr = vertcat(freqCorr{:});

for ii=1:size(totRelInds, 1)
    protocolStruct.stim(ii).matCell = totMats{ii};
    protocolStruct.stim(ii).relInds = totRelInds(ii,:);
    protocolStruct.stim(ii).freqCorr = totFreqCorr(ii);
end


fprintf('Protocol structure created \n')




end
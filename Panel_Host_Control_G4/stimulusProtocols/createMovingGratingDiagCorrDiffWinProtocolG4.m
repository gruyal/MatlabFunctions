function protocolStruct = createMovingGratingDiagCorrDiffWinProtocolG4(inputStruct)

% function createMovingGratingDiagCorrDiffWinProtocolG4(inputStruct)
%
% This function uses the inputStruct and createProtocol function to
% generate a stationary grating that will appear in all the phases. It has certain assumptions and therefore requires
% less inputs.
% If function is called with no inputStruct, it prompt the user for
% relevant input and allows to change the default
%
% This function generates a moving grating in both directions of the given orientation
% It is similar to movingGrating only instead of changing phase it changes
% window size
%
% ASSUMPTIONS
% Randomize -       everything is randomized (so that if several gratings appear
%                   they appear in all positions and in all oreintations.
% Masks positions - [1,1] Grid is assumed.
% gratingFuncHand - Function uses the generateGrstingFrameByInds to make sure
%                   even in diagonal the coverage is complete
% maskType -        Due to the way the bar is generated, mask type is
%                   always rectangle.
% grtMaskInt -      one to one grating and mask. Set to 1. Since mask need to change orientation also
%
%   NOTE!!! to correct small problems with diagonal line orientation is
%   implemented here and not in createProtocol
%
% INPUT
% Defaults are dilimited with {} and are optional inputs.
% Only 2 obligatory fields are width and gridCenter
%
% inputStruct -     Should have the following fields
% .f/sBarV -        1XN Vector. normalized luminance value for first and second bars in the grating. { 1 / 0  }
%                   Note!! if length for fBar and sBar should be identical
% .winSize -        1XM vector. size of grating window in pixels. Since the
%                   protocol uses different sized windows, this size
%                   applies for both the width and the height of the window
%                   (will always be square - though rectangle is used to
%                   genrate the window)
% .grtHWidth -      1XW vector(in pixels). Width of a bar in the grating,
%                   effectively half the cycle { 4 }.
%                   For grating with width w, there are 2*w phases
% .orientation -    single number (0-3). Since the bar isn't moving 4-7 are redundant.
%                   Orientations for the gratings. Applied on all inputs.
% .numCyc -         Single number. Number of cycles to present (a cycle is
%                   all the phases of a grating) { 3 }
% .stepDur -        1XN vector. Duration in seconds for each phase of the grating
% .gsLevel -        gray scale level ( fixed at 4 )
% .gratingMidVal -  value of the rest of the window (0.49 - bkgd level)
% .gridCenter -     1X2 vector specifying the center of the grid in X and Y
%                   (sptial coordinates in pixels <for an 8X4 arena its 96X32). If one dimension of
%                   grid is even, grid will be presented around center but
%                   will not have a position in the actual center.
% .intFrames -      number of empty intervening frames. If not given half a
%                   second worth (based on generalFrequency)
% .repeats -        scalar. number of times the whole protocol repeats (passed into createProtocol) {3}
% .generalFrequency-Frequency with which frames from the generated protocol
%                   will be dumped (passed on to runDumpProtocol) in position function units
%                   (frames per second on the controller). fixed at 500 for gsLevel 4
% .freqCorrFlag -   Also passed on to runDumpProtocol. Logical flag to
%                   indicate whether different stimuli should be run with temporal frequency
%                   correction { 1 }.
%
% OUTPUT
% protocolStruct with all the required fields from createProtocl.
%
%
%                   NOTE! masks and gratings need not be of the same length
%                   NOTE! grid mask positions will use only the first
%                   maskRadius value to interpert overlap.

%% GENERAL AND DEFAULT PARAMETERS

baseSiz = 445; % size of single frame or mask
arenaSize = [192,48];
gratingFuncHand = @generateGratingFrameByInds;

default.fBarV = 1;
default.sBarV = 0;
default.winSize = [9, 15, 21];
default.grtHWidth = 4;
default.numCyc = 5;  % 3 was too short
default.gridCenter = 'UI';
default.gratingMidVal = 0.49;
default.orientations = 'UI';
default.stepDur = [0.04, 0.16];
default.intFrames = nan;
default.repeats = 3;
default.randomize = 1;

fixed.gsLevel = 4;
fixed.gridSize = [1,1];
fixed.gridOverlap = 0;
fixed.generalFrequency = 500;
fixed.maskType = {'rectangle'};
fixed.freqCorrFlag = 0;
fixed.grtMaskInt = 1;

% combining default and input structures
if nargin == 0
    default = modifyDefaultStruct(default);
else
    default = modifyDefaultStruct(default, inputStruct);
end


 %% ORIENTATIONS

 ort = default.orientations;
 assert(length(ort) == 1, 'Orientation should be single number')
 assert(ismember(ort, 0:3), 'Orientation values should be between 0 and 3')

 newOrt = ort; %(ismember(ort, 0:3));


 % orientation is implemented internally
 protocolStruct.orientations = 0;




  %% MASK (masks created with grating)

 maskHS = floor(default.winSize/2);
 assert(isvector(maskHS), 'winHeight should be a 1XM vector');
 assert(min(maskHS) > 0, 'winHeight minimum should be a positive number')

 maskT = fixed.maskType;

%% GRATING PARAMETERS

% needed to determine number of frames to appear
protocolStruct.generalFrequency = fixed.generalFrequency;

stepLen = sort(default.stepDur);
assert(isvector(stepLen), 'stepDur should be one a vector')

stepFrames = unique(round(stepLen * fixed.generalFrequency));
if length(stepFrames) < length(stepLen)
    warning('%d step durations omitted since were the same after rounding', length(stepLen) - length(stepFrames))
end

assert(min(stepFrames) > 0, 'stimulus can not be presented for such a short duration. Minimal duration is 20ms')

gsLev = fixed.gsLevel;

bkgdVal = default.gratingMidVal;
assert(length(bkgdVal) == 1, 'gratingMidVal should be a single number');
assert(bkgdVal >=0 && bkgdVal <= 1, 'gratingMidVal should be between 0 and 1')

stimFB = default.fBarV;
assert(isvector(stimFB), 'fBarV should be a 1XN vector')
assert(min(stimFB) >=0 && max(stimFB) <= 1, 'fBarV should be between 0 and 1');

stimSB = default.sBarV;
assert(isvector(stimSB), 'sBarV should be a 1XN vector')
assert(min(stimSB) >=0 && max(stimSB) <= 1, 'sBarV should be between 0 and 1');

assert(length(stimSB) == length(stimFB), 'fBarV and sBarV should have the same length')


numCyc = default.numCyc;
assert(length(numCyc) == 1 , 'numCyc should be a single number')
assert(numCyc == round(numCyc), 'numCyc should be integer');
assert(numCyc >= 1, 'numCyc should be positive number');

barW = default.grtHWidth;
assert(isvector(barW), 'barWidth should be a vector');
assert(all(bsxfun(@eq, barW, round(barW))), 'barWidth should use only integers');

count = 0;
gratingArray = [];

for hh=1:length(maskHS)

    relRegR = maskHS(hh);
    relDiagR = round((2*relRegR+1)/sqrt(2)); % was +1 (with -1 overlap with non-rotated square is too small);

    radVec = [relRegR, relDiagR];

    ortPosInd = rem(newOrt,2)+1;
    relRad = radVec(ortPosInd);

    sqDim = 2*maskHS(hh)+1;


    for ss=1:length(stimFB)

        for kk=1:length(stepFrames)

            for ww=1:length(barW)

                relPhase = 1:2*barW(ww);
                phaseMat = zeros(2, length(relPhase));
                phaseMat(1,:) = relPhase;
                phaseMat(2,:) = circshift(fliplr(relPhase),[0, 1]);
%                 phaseMat(3,:) = circshift(relPhase, [0, barW(ww)]);
%                 phaseMat(4,:) = circshift(fliplr(relPhase), [0, barW(ww)+1]);

                for dd=1:size(phaseMat,1)

                    basePhase = phaseMat(dd, :);
%                     basePhase = reshape(repmat(preBasePhase, stepFrames(kk), 1), 1, []);
                    corrPhase = repmat(basePhase, 1, numCyc);

                    count = count+1;
                    gtStruct(count).wid = barW(ww);
                    gtStruct(count).ori = newOrt;
                    gtStruct(count).fVal = stimFB(ss);
                    gtStruct(count).sVal = stimSB(ss);
                    gtStruct(count).sqDim = sqDim; % span is a bit unique in its claculation here
                    gtStruct(count).phase = corrPhase;
                    gtStruct(count).gsLevel = gsLev;
                    gtStruct(count).bkgdVal = bkgdVal;
                    gtStruct(count).matSize = baseSiz;
                    gtStruct(count).stepFrames = stepFrames(kk); 

                    maskSt(count).type = maskT{1};
                    maskSt(count).radius = [relRad, maskHS(hh)];
                    maskSt(count).ori = newOrt;

                    direc = basePhase(3) - basePhase(2); % since there is sometimes a jump between first and second
                    stPhase = basePhase(1);
                    if stPhase > 1
                        stPhase = -1; % to make it comparable between diff widths
                    end

                    gratingArray = vertcat(gratingArray, ...
                                           [count, stimFB(ss), stimSB(ss), barW(ww), newOrt, 2*maskHS(hh)+1, ...
                                            direc, stepFrames(kk)/fixed.generalFrequency]);

                end
            end
        end
    end
end

 tabVarNames =  {'index','FBval', 'SBVal', 'width', 'orient', 'WinSize', 'direction', 'stimDur'};

 protocolStruct.gratingTable = array2table(gratingArray, 'variablenames', tabVarNames);
 protocolStruct.gratingStruct = gtStruct;
 protocolStruct.masksStruct = maskSt;
 protocolStruct.relGtStName = 'phase';
 
 %% GRID

% since the grid is just one position

gdCen = default.gridCenter;
assert(isvector(gdCen), 'gridCenter should be a 1X2 vector')
assert(length(gdCen)==2, 'gridCenter should be a 1X2 vector')

gridSt.gridSize = fixed.gridSize;
fixed.gridOverlap;
gridSt.spacing = [maskHS(1), maskHS(1)];
gridSt.startPos = gdCen;

maskPos = makeGrid(gridSt);

protocolStruct.maskPositions = maskPos;


 %% Misc parameters

 protocolStruct.freqCorrFlag = fixed.freqCorrFlag;

 protocolStruct.funcHand = gratingFuncHand;
 protocolStruct.interleave = fixed.grtMaskInt;

 intF = default.intFrames;
 if isnan(intF)
     protocolStruct.intFrames = floor(fixed.generalFrequency/2);
 else % if user gave a number
    assert(intF >= 0, 'intFrames should be a non-negative number')
    protocolStruct.intFrames = intF;
 end

 protocolStruct.repeats = default.repeats;

 protocolStruct.randomize.gratingSeq = default.randomize;
 protocolStruct.randomize.masks = default.randomize;
 protocolStruct.randomize.orientations = default.randomize;
 protocolStruct.randomize.maskPositions = default.randomize;

 %% Creating protocl


 protocolStruct = createProtocolG4(protocolStruct);

 protocolStruct.inputParams = default;


end

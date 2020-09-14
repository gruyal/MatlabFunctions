function protocolStruct = createMovingBarDiagCorrProtocol(inputStruct)

% function createMovingBarDiagCorrProtocol(inputStruct)
%
% This function uses the inputStruct and createProtocol function to
% generate single bar moving through the window . It has certain assumptions and therefore requires
% less inputs.
% If function is called with no inputStruct, it prompt the user for
% relevant input and allows to change the default
%
% ASSUMPTIONS
% Randomize -       everything is randomized (so that if several gratings appear
%                   they appear in all positions and in all oreintations.
% Masks positions - [1,1] Grid is assumed.
% gratingFuncHand - Function uses the generateBarFrameByInd to make sure
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
% .barWid -         1XW vector of width of bar in pixels
% .stimBar -        1XN vector (0-1) normalized luminance value.
% .barHeight -      1XM vector. height of bar in pixels.
% .barSpan -        1XS vector. span along which bar will move (in
%                   pixels). Will be converted into the width of the
%                   rectangular mask
% .orientation -    Vector (0-3). Since the bar isn't moving 4-7 are redundant.
%                   Orientations for the gratings. Applied on all inputs.
% .stepDur -        duration in seconds in which the bar will appear
% .gsLevel -        gray scale level ( fixed at 4 )
% .gratingMidVal -  value of the rest of the window (0.49 - bkgd level)
% .maskPositions -  User can specify these directly as an NX2 matrix, or
%                   use the other parameters to generate them (if this is
%                   given other parameters are disregarded
% .gridSize  -      1X2 vector specifying size of grid in X and Y (spatial
%                   coordinates). { [2,2] }
% .gridPosOverlap - Overlap between different positions on mask position grid.
%                   Units s are normalized maskSizes so that 0 means no overlap and no gap,
%                   1 means complete overlap (meaningless), and -1 means a gap of one mask
%                   between positions. { 0 }
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
gratingFuncHand = @generateBarFrameByInds;

default.stimBar = 1;
default.barWid = [1,2];
default.barHeight = 9;
default.barSpan = 9;
default.gridCenter = 'UI';
default.gratingMidVal = 0.49;
default.orientations = 0:7;
default.stepDur = 0.08;
default.intFrames = nan;
default.repeats = 3;
default.randomize = 1;


fixed.gsLevel = 4;
fixed.generalFrequency = 500;
fixed.maskType = {'rectangle'};
fixed.freqCorrFlag = 0;
fixed.grtMaskInt = 1;
fixed.gridSize = [1,1];
fixed.gridOverlap = 0;

% combining default and input structures
if nargin == 0
    default = modifyDefaultStruct(default);
else
    default = modifyDefaultStruct(default, inputStruct);
end


 %% ORIENTATIONS

 ort = default.orientations;
 assert(isvector(ort), 'Orientation should be 1XM vector')
 assert(prod(ismember(ort, 0:7)) == 1, 'Orientation values should be between 0 and 7')

 newOrt = ort; %(ismember(ort, 0:3));


 % orientation is implemented internally
 protocolStruct.orientations = 0;




  %% MASK (masks created with grating)

 maskHW = floor(default.barSpan/2); % rectangle mask input is half width
 assert(isvector(maskHW), 'barSpan should be a 1XS vector')
 assert(min(maskHW) > 0, 'barSpan should be positive')


 minMaskR = min(maskHW);

 maskHH = floor(default.barHeight/2);
 assert(isvector(maskHH), 'barHeight should be a 1XM vector');
 assert(min(maskHH) > 0, 'barHeight minimum should be a positive number')

 maskT = fixed.maskType;

%  maskCount = 0;
%  for ii=1:length(maskHH)
%      maskCount = maskCount+1;
%      maskSt(maskCount).type = maskT{1};
%      maskSt(maskCount).radius = [maskHW, maskHH(ii)];
%  end
%
%
%  numMasks = length(maskSt);
%
%  protocolStruct.masksStruct = maskSt;

%% GRATING PARAMETERS

% needed to determine number of frames to appear
protocolStruct.generalFrequency = fixed.generalFrequency;

stepLen = default.stepDur;
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

stimB = default.stimBar;
assert(min(stimB) >=0 && max(stimB) <= 1, 'stimBar should be between 0 and 1');
assert(isvector(stimB), 'stimBar should be a 1XN vector')

barW = default.barWid;
assert(isvector(barW), 'barW should be a 1XW vector')
assert(min(barW) > 0, 'barW should be a positive number')


count = 0;
gratingArray = [];

for vv=1:length(stimB)

    for ww=1:length(barW)

        for hh=1:length(maskHH)

            for sp = 1:length(maskHW)

                relRegR = maskHW(sp);
                relDiagR = round((2*relRegR+1)/sqrt(2));% was +1 (with -1 overlap with non-rotated square is too small);
                radCell = {relRegR, relDiagR};

                for kk=1:length(stepFrames)

                    for oo=1:length(newOrt)

                        count = count+1;

                        ortPosInd = rem(newOrt(oo),2) +1;
                        relRad = radCell{ortPosInd};

                        relPos = -relRad:relRad+barW(ww)-1;
%                         corrPos = reshape(repmat(relPos, stepFrames(kk), 1), 1, []);

                        gtStruct(count).wid = barW(ww);
                        gtStruct(count).ori = newOrt(oo);
                        gtStruct(count).val = stimB(vv);
                        %gtStruct(count).sqDim = max(2*maskHW(sp)+1, 2*maskHH(hh)+1); % generateBarFrameByInds corrects for diagonal internally
                        gtStruct(count).sqDim = 2*maskHW(sp)+1; % since when using divideTotSquareToCols height is not taken into account

                        gtStruct(count).pos = relPos;
                        gtStruct(count).stepFrames = stepFrames(kk);
                        gtStruct(count).gsLevel = gsLev;
                        gtStruct(count).bkgdVal = bkgdVal;
                        gtStruct(count).matSize = baseSiz;

                        maskSt(count).type = maskT{1};
                        maskSt(count).radius = [relRad, maskHH(hh)];
                        maskSt(count).ori = newOrt(oo);

                        gratingArray = vertcat(gratingArray, ...
                                                    [count, stimB(vv), barW(ww), 2*maskHW(sp)+1, ...
                                                    stepFrames(kk)/fixed.generalFrequency, 2*maskHH(hh)+1, newOrt(oo)]);

                    end
                end
            end
        end
    end
end


 tabVarNames =  {'index', 'value', 'width', 'span', 'stepDur', 'height', 'orient'};

 protocolStruct.gratingTable = array2table(gratingArray, 'variablenames', tabVarNames);
 protocolStruct.gratingStruct = gtStruct;
 protocolStruct.masksStruct = maskSt;

 %% GRID

 gridSt.gridSize = fixed.gridSize;
 ovlp = fixed.gridOverlap;

 maskS = 2*minMaskR+1;
 space = maskS - maskS*ovlp;
 gridSt.spacing = [space, space];

 gdCen = default.gridCenter;
 assert(isvector(gdCen), 'gridCenter should be a 1X2 vector')
 assert(length(gdCen)==2, 'gridCenter should be a 1X2 vector')

 stCrds = gdCen - space*(gridSt.gridSize-1)/2;
 if stCrds(1) < 1
    warning('Grid start position in X is out of range - changed to 1')
    stCrds(1) = 1;
 end
 if stCrds(2) < 1
    warning('Grid start position in Y is out of range - changed to 1')
    stCrds(2) = 1;
 end

 gridSt.startPos = stCrds;

 maskPos = makeGrid(gridSt);

 protocolStruct.maskPositions = maskPos;



 %% Misc parameters

 protocolStruct.freqCorrFlag = fixed.freqCorrFlag;

 protocolStruct.funcHand = gratingFuncHand;
 protocolStruct.interleave = fixed.grtMaskInt;

 intF = default.intFrames;
 if isnan(intF)
     protocolStruct.intFrames = floor(fixed.generalFrequency/2.5);
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

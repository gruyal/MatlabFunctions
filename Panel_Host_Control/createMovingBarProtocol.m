function protocolStruct = createMovingBarProtocol(inputStruct)

% function createFlickerBarProtocol(inputStruct)
%
% This function uses the inputStruct and createProtocol function to
% generate a long narrow bar that moves across the positions of the window. 
% The stimulus is designed to map the center of a receptive field relatively quickly. 
% It has certain assumptions and therefore requires less inputs. 
% If function is called with no inputStruct, it prompt the user for
% relevant input and allows to change the default
% 
% ASSUMPTIONS
% Since this is a grating certain parameters are assumed to be symmetrical.
% 
% position -        Bar will move across in all window positions
% Contrast -        each bar is uniform (no on/off gradients), and both are same
%                   distance from mid level GS (background). 
% Randomize -       everything is randomized (so that if several gratings appear
%                   they appear in all positions and in all oreintations. 
% grtMaskInt -      since grating width is completely dependent on mask
%                   size, they will be read one to one. {interleave 1}
% gratingFuncHand - Function uses the generateGratingFrame
%
%
% INPUT 
% Defaults are dilimited with {} and are optional inputs. 
% Only 2 obligatory fields are maskRadius and gridCenter
%
% inputStruct -     Should have the following fields
% .barWidth  -      1XN vector. Width of the moving bar. { 2 }
% .barCont -        at most a 2 element vector. contrast of bar moving. For now only 0 or 1 (can be
%                   adjusted by changing contract. {[0,1]}
% .gsLevel -        gray scale level for grating frames
% .contrast -       Difference between bright and dark bars relative to
%                   mid. For now length should be 1.
% .maskPositions -  User can specify these directly as an NX2 matrix, or
%                   use the other parameters to generate them (if this is
%                   given other parameters are disregarded
% .maskType -       'circle' or {'square'}. Can also be a 1XM cell array
% .maskRadius -     1XP vector in pixels. mask is actually 2XmaskRadius+1. {2*width - closest value from the available ones} { 4 }   
% .maskInt -        logical. If TRUE function will generate all the combinations
%                   between type and radius. { 17 } 
% .gridSize  -      1X2 vector specifying size of grid in X and Y (spatial
%                   coordinates). { [1,1] }
% .gridPosOverlap - Overlap between different positions on mask position grid. 
%                   Units s are normalized maskSizes so that 0 means no overlap and no gap, 
%                   1 means complete overlap (meaningless), and -1 means a gap of one mask
%                   between positions. { 0 }
% .orientations -   bar movement direction { 0:2:6 - 4 cardinal directions}
% .gridCenter -     1X2 vector specifying the center of the grid in X and Y
%                   (sptial coordinates in pixels <for an 8X4 arena its 96X32). If one dimension of
%                   grid is even, grid will be presented around center but
%                   will not have a position in the actual center. {[48, 16]}
% .intFrames -      number of empty intervening frames. If not given fifth a
%                   second worth (based on generalFrequency)
% .repeats -        scalar. number of times the whole protocol repeats (passed into createProtocol) {3} 
% .generalFrequency-Frequency with which frames from the generated protocol
%                   will be dumped (passed on to runDumpProtocol) in Hz.
%                   { 15 }
% .freqCorrFlag -   Also passed on to runDumpProtocol. Logical flag to
%                   indicate whether different stimuli should be run with temporal frequency
%                   correction { 0 }.  
%
% OUTPUT 
% protocolStruct with all the required fields from createProtocl. 
%
%
%                   NOTE! masks and gratings need not be of the same length
%                   NOTE! grid mask positions will use only the first
%                   maskRadius value to interpert overlap. 
%                   Note! since mask is bigger than arena's height, first
%                   and last position of horizontal bars are not displayed 

%% GENERAL AND DEFAULT PARAMETERS

baseSiz = 225; % size of single frame or mask
arenaSize = [96, 32];
gratingFuncHand = @generateGratingFrame;

default.gridCenter = [48,16];
default.generalFrequency = 15;
default.barWidth = 2; 
default.barCont = [0,1];
default.contrast = 1;
default.orientations = 0:2:6;
default.gsLevel = 3;
default.maskInt = 1;
default.gridOverlap = 0;
default.grtMaskInt = 1;
default.gratingMidVal = 0.49;
default.intFrames = nan;
default.repeats = 3;
default.randomize = 1;
default.freqCorrFlag = 0;


% fixed parameters for this protocol that are not user modified
fixed.maskRadius = 17;
fixed.maskType = {'square'};
fixed.gridSize = [1,1];

% combining default and input structures
if nargin ==0
    default = modifyDefaultStruct(default);
else
    default = modifyDefaultStruct(default, inputStruct);
end


 %% ORIENTATIONS
 
 ort = default.orientations;
 assert(isvector(ort), 'Orientation should be 1XM vector')
 assert(prod(ismember(ort, 0:7)) == 1, 'Orientation values should be between 0 and 7')
 if sum(ismember(ort, 1:2:7) > 0)
     warning('odd orientations will give a diagonal window')
 end
 protocolStruct.orientations = ort;
 
 %% MASK
 
 maskT = fixed.maskType;
 maskR = fixed.maskRadius;
 
 maskSt(1).type = maskT{1};
 maskSt(1).radius = maskR(1);
 
 % might change after startBar is read in 
 numMasks = length(maskSt);

 %% GRATING PARAMETERS
 
 barW = default.barWidth;
 assert(isvector(barW), 'barWidth should be a vector')
 
 barC = default.barCont;
 assert(isvector(barC) && length(barC) <= 2, 'barCont should be a vector of at most length 2')
 assert(prod(ismember(barC, [0,1])) == 1, 'barCont can only be 0, 1, or both')
 
 cont = default.contrast;
 assert(isvector(cont) && length(cont) == 1, 'Contrast should be vector of length 1')
 assert(cont > 0 && cont <= 1, 'Contrast should be  0 < cont <= 1');
 
 bkgdVal = default.gratingMidVal;
 assert(length(bkgdVal) == 1, 'gratingMidVal should be a single number');
 assert(bkgdVal >=0 && bkgdVal <= 1, 'gratingMidVal should be between 0 and 1')
 
 onVal = bkgdVal + cont/2; % 0.49 is for the middle value to be rounded down (in GS3 it is 3 and not 4)
 offVal = bkgdVal - cont/2.041; % so that it wont go negative
 
 
 gsLev = default.gsLevel;
 assert(ismember(gsLev, 1:4), 'gsLevel should be an integer between 1 and 4')
 
 count = 0;
 relMaskR = maskSt(1).radius;
 
for ii=1:length(barW)
       
    for jj=1:length(barC)
        
        count = count+1;
        if barC(jj) == 0
            gtStruct(count).valsONSt = bkgdVal;
            gtStruct(count).valsONEnd = bkgdVal;
            gtStruct(count).valsOFFSt = offVal;
            gtStruct(count).valsOFFEnd = offVal;
            gtStruct(count).widthON = 2*relMaskR+1; % to make sure bar gradually appears and disappears
            gtStruct(count).widthOFF = barW(ii);
            gtStruct(count).barAtPos = 0;
        elseif barC(jj) == 1
            gtStruct(count).valsONSt = onVal;
            gtStruct(count).valsONEnd = onVal;
            gtStruct(count).valsOFFSt = bkgdVal;
            gtStruct(count).valsOFFEnd = bkgdVal;
            gtStruct(count).widthOFF = 2*relMaskR+1;
            gtStruct(count).widthON = barW(ii);
            gtStruct(count).barAtPos = 1;
        end
        
        gtStruct(count).position = -(relMaskR+barW(ii)):relMaskR;
        gtStruct(count).gsLevel = gsLev; 
        
        newMaskSt(count) = maskSt; % since length(maskSt) == 1  

    end
    
 end
 
 protocolStruct.gratingStruct = gtStruct;
 protocolStruct.masksStruct = newMaskSt;
 
 %% GRID 
 
 % if maskPositions exist then the rest of the parameters are unnecessary
 if nargin  ==1 && isfield(inputStruct, 'maskPositions')
     maskPos = inputStruct.maskPositions;
     assert(size(maskPos, 2) == 2, 'mask positions should be an NX2 matrix')
     assert(max(maskPos(:,1)) <= arenaSize(1), 'mask position X values should not exceed %d', arenaSize(1))
     assert(max(maskPos(:,2)) <= arenaSize(2), 'mask position Y values should not exceed %d', arenaSize(2))
     assert(min(maskPos(:)) > 0, 'mask position values should be positive')
     protocolStruct.maskPositions = maskPos;
     
 else %if maskPositions isn't given use default paramters (or grid input)
     
     gridSt.gridSize = fixed.gridSize;
     
     ovlp = default.gridOverlap;
     assert(isscalar(ovlp), 'Overlap should be a single number')
     assert(ovlp < 1, 'Overlap value above 1 are not accepted')
     assert(ovlp ~=1, 'What are you stupid? overlap 1 means no grid')
 
     maskS = 2*maskR(1)+1;
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
    
 end
 
 %% Misc parameters
 
 protocolStruct.freqCorrFlag = default.freqCorrFlag;
 protocolStruct.generalFrequency = default.generalFrequency;
 
 protocolStruct.funcHand = gratingFuncHand;
 protocolStruct.interleave = default.grtMaskInt;
 
intF = default.intFrames;
 if isnan(intF)
     protocolStruct.intFrames = floor(default.generalFrequency/5);
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
 
 protocolStruct = createProtocol(protocolStruct);
 
 protocolStruct.inputParams = default;
 
end











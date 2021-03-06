function makePatandPosFuncForBExp_barExpand(stimStruct)

% function makePatternandPosFuncForBehavior(stimStruct)
%
% This function generates patterns and position functions for bheavioral
% experiments (oscilating an object at a specified location with specified
% amp). It will generate the pattern and posFunc file in the working
% directory naming them based on parameters in expStruct. 
% The first dimension in the inputs (center, arenaSize etc.) refers to the
% vertical dimension in the physical arena from top to bottom. The second
% refers to rhe horizontal dimension from left to right
%
% INPUT
%
% stimStruct -           structure that should have certain fields 
%       .barSt:         structure with the necessary fields for placeBarInArena
%       .halfAmp:       amplitude of movement in pixels in the horizontal dimension
%       .center:        (optional) 1X2 vector for center of motion. If not given
%                       center of arena is used
%       .background:    (optional) string describing the required backrround for the bar. 
%                       Options are: {'uniform'} <mid level GS>, 
%                       'randomS' <random dots static>, 
%                       'randomC' <random dot changing each frame> 
%       .bkgdProp:      (optional) 1X3 vector. if one of the random options is used bkgd prop will
%                       be used to generate the background (the same way that proportion is
%                       used in placeBarInArena. Default is [1/3, 1/3, 1/3]
%       .barSize:       outside of the barSt, this parameter should be an Nx2 matrix
%                       with each row specifying x and y of bar size (would be used as .size in
%                       placeBarInArena. 
%       .barStable      (optional) logical. if {TRUE} pixels will not change relative
%                       position with bar movement. 
%       .expand         (optional) string. determines whether the expanding bar will
%                       preserve the {'number'} or 'proportion' of pixels. When number is
%                       chosen, the added pixels will be of mid-Value (see placeBarinArena)
%     
%
%   NOTE: First element should be the smallest and number of pixels that
%   are ON/OFF would be maintained for the rest of the bar sizes
%   
% OUPUT
% pattern structure with position function at specified freqeuncies
% generated for it and embedded as one of the fields. pattern fields are:
%
% .Pats, .x_num, .y_num, .num_panels, .gs_val, .BitMapIndex, .data 
% all like in the previous use of the pattern structure with the original controller 
% .stimStruct:          structure that was used to generate the pattern
% .stimMap:             container.Map object that specifies characteristics of the
%                       pattern (for easier access)
% .freqVecMap:          map object with keys being frequnecy and values
%                       actual position functions
% 
% % Perhaps add background range in the future


arenaSize = [32, 96]; %in pixels
gsLev = 3;
midVal = floor((2^gsLev-1)/2);


% checking input
if ~isfield(stimStruct, 'barSt')
    stimStruct.barSt = generateDefaultBarSt(arenaSize);
end
assert(isfield(stimStruct, 'halfAmp'), 'expStruct is missing halfAmp field')
assert(length(stimStruct.halfAmp) == 1, 'halfAmp should be a single number')
halfAmp = stimStruct.halfAmp;

if isfield(stimStruct, 'center')
    assert(stimStruct.center(2) + stimStruct.halfAmp <= arenaSize(2), 'oscilation is out of range')
    assert(stimStruct.center(2) - stimStruct.halfAmp >= 1, 'oscilation is out of range')
    center = stimStruct.center;
else
    center = [ceil(arenaSize(1)*2/3), floor(arenaSize(2)/2)];
end

if isfield(stimStruct, 'background')
    bkgdStr = lower(stimStruct.background);
    assert(ismember(bkgdStr, {'uniform', 'checker', 'randoms', 'randomc'}), ...
           'background can be only one of the following: uniform, checker, randoms, randomc')
    bkgd = bkgdStr;   
else
    bkgd = 'uniform';
end

assert(isfield(stimStruct, 'barSize'), 'barSize is missing from expStruct')
assert(ismatrix(stimStruct.barSize), 'barSize should be an Nx2 matrix')
assert(size(stimStruct.barSize, 2) == 2, 'barSize should be an Nx2 matrix')

if isfield(stimStruct, 'barStable')
    barS = stimStruct.barStable;
else
    barS = 1;
end

if isfield(stimStruct, 'expand')
    exStr = lower(stimStruct.expand);
    assert(ismember(exStr, {'number', 'proportion'}), 'expand should be either number or proportion')
    exVal = exStr;
else
    exVal = 'number';
end


relRange = -halfAmp:halfAmp;
bkgdChange = 0;


% creating background image
if      strcmp(bkgd, 'uniform')
    bkgdIm = ones(arenaSize)*midVal;
elseif  strcmpi(bkgd, 'checker')
    bkgdIm = generateCheckboard(arenaSize)*(2^gsLev-1);
elseif  strcmp(bkgd, 'randoms')
    if isfield(stimStruct, 'bkgdProp')
        prop = stimStruct.bkgdProp;
    else
        prop = [1, 1, 1];
    end
    bkgdIm = generateRandomBkgd(arenaSize, prop);
    bkgdIm(bkgdIm == 1) = midVal;
    bkgdIm(bkgdIm == 2) = 2^gsLev-1;
elseif strcmp(bkgd, 'randomc')
    bkgdChange = 1;
    if isfield(stimStruct, 'bkgdProp')
        prop = stimStruct.bkgdProp;
    else
        prop = [1, 1, 1];
    end
    bkgdIm = ones([arenaSize, length(relRange)]);
    for ii=1:size(bkgdIm, 3)
        bkgdIm(:,:,ii) = generateRandomBkgd(arenaSize, prop);
    end
    bkgdIm(bkgdIm == 1) = midVal;
    bkgdIm(bkgdIm == 2) = 2^gsLev-1;
end
    
% adjusting bar so that number of pixels does not change with size

minObjSize = min(prod(stimStruct.barSize,2));

if strcmp(exVal, 'number')
    relProp = stimStruct.barSt.proportion/sum(stimStruct.barSt.proportion);
    numMinPix = relProp(1) * minObjSize;
    numMaxPix = relProp(3) * minObjSize;
end

tempBarSt = stimStruct.barSt;
basicStim = ones([arenaSize, length(relRange)+1, size(stimStruct.barSize, 1)+1]) * midVal; % leaves an empty frame in the beginning 

for ii = 1:size(stimStruct.barSize,1)
    tempBarSt.size = stimStruct.barSize(ii, :);
    if strcmp(exVal, 'number')
        tempBarSt.proportion = [numMinPix, prod(tempBarSt.size) - numMinPix - numMaxPix, numMaxPix];
    elseif strcmp(exVal, 'proportion')
        tempBarSt.proportion = stimStruct.barSt.proportion;
    end
    rngSeed = now;
    
    for jj=1:length(relRange)
        tempBarSt.center = center + [0, relRange(jj)];
        if barS % stabalizes the bar even if it is randomly generated
            tempBarSt.seed = rngSeed;
        end
        
        if bkgdChange
            basicStim(:, :, jj+1, ii+1) = placeBarInArena(bkgdIm(:,:,jj), tempBarSt);
        else
            basicStim(:, :, jj+1, ii+1) = placeBarInArena(bkgdIm, tempBarSt);
        end
    end
end


%% Creating pattern structure and posFunc files


pattern.Pats        = basicStim;
pattern.x_num       = size(basicStim,3);
pattern.y_num       = size(basicStim,4);
pattern.num_panels  = 48;
pattern.gs_val      = gsLev;     % 8 levels of intensity (0-7)

panel_id_map =                  [ 4 11  7  3 10  6  2  9  5  1 12  8 ; 
                                 16 23 19 15 22 18 14 21 17 13 24 20 ;
                                 28 35 31 27 34 30 26 33 29 25 36 32 ;
                                 40 47 43 39 46 42 38 45 41 37 48 44 ];%for fly treadmill id
pattern.Panel_map = panel_id_map;
pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = Make_pattern_vector(pattern);
pattern.stimStruct = stimStruct; % to have a copy of the structure with the saved pattern

stimMapObj = getStimStructMap(stimStruct);
stimMOLab = cellfun(@(x) [num2str(x), '_'], values(stimMapObj), 'uniformoutput', 0); % to be used as string for pattern name
stimMOLab = [stimMOLab{:}];
stimMOLab = stimMOLab(1:end-1); % gets rid of last underscore 

pattern.stimMap = stimMapObj;


%% save pattern file
bExpTitle = 'barExpand_';

fName = ['PatNPosFunc_', bExpTitle, stimMOLab, '.mat'];


% Makeing position functions

freqVec = 4:4:20;
xpos = pattern.x_num;
% funcName = ['position_function_X_', bExpTitle, stimMOLab, '_freq_']; 
posFuncMapObj = containers.Map('keyType', 'double', 'valueType', 'any');

for ii=1:length(freqVec)
    
    tempFunc = makeSinPosFunc(xpos, freqVec(ii));
    posFuncMapObj(freqVec(ii)) = tempFunc;
    
end

pattern.freqVecMap = posFuncMapObj;

% overwrite protection
fn = dir('Pattern*.mat');
fnCell = {fn.name};

if ~isempty(fnCell)
    if ismember(fName, fnCell)
        butVal = questdlg('File already exists. Overwrite?', 'File Error', 'Yes');
        
        switch butVal
            case 'Yes'
                save(fName, 'pattern')
            case 'No'
                disp('file saved with time stamp')
                save([fName(1:end-4), '_', datestr(now, 'MMSS'), '.mat'], 'pattern')
            case 'Cancel'
                disp('no file saved!')
                return
        end
    else
       save(fName, 'pattern') 
        
    end
else
    save(fName, 'pattern')
end





end

%%

function randBkgd = generateRandomBkgd(arenaSize, prop)
% generates background image based on the given proportions. 0 for min
% values, 1 for mid val, and 2 for max vals

baseProp = floor(prop/sum(prop) * prod(arenaSize));    
randBkgd = ones(arenaSize);
minInd = randsample(1:prod(arenaSize), baseProp(1));
indsLeft = setdiff(1:prod(arenaSize), minInd);
maxInd = randsample(indsLeft, baseProp(3));

randBkgd(minInd) = 0;
randBkgd(maxInd) = 2;

end

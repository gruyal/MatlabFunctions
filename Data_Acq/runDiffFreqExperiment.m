function expResults = runDiffFreqExperiment(direct, iniPos, protInds, freqVec, random)

% function data = runPanelExperiment(direct)
% This function runs the protocols generated by 'createExpStruct' and
% 'createProtocolsScript' and collects the gennerated data. 
% The function runs a session for each seperate protocol and saves after each one 
%
% INPUT
% direct -      relevant 'on_SDcard' directory in which expStruct is save and to which 
%               the generated data would be saved. expStruct is generated as
%               stated above, and should contain pattern, posFunc, and protocol
%               fields. It is also recommended to use 'checkExpStruct' to verify
%               that the protocols are referring to the right files. 
% iniPos -    offset to apply for protocols labeled with a [-1 1] initial position. 
%               in this version it is a single value for the x position
%               alone 
% protInds -    Vector of relevant protocol indices. if not given
%               the function will use all protocols defined as 'short' in 'type' field.
% freqVec -     Vector of different freqeuncies with which to run to
%               relevant protocols
% random -      logical (optional). Flag for whether the protocols should be 
%               presented in a random order or not. Default is random
%
% OUTPUT 
% shortExpInfo  -   saves several variable that are related to the experiment
% expResults -      data aquired for each protocol is an N X 5 matrix with N being 
%                   the number of samples, and 5 channels which are: 
%                   (1) membrane Current (2) membrane Voltage (3)X position (4) Y
%                   position (5)photodiode data
% expStruct -       The actual expStructure used to run the experiment
% expProt -         experimental protocols ordered by presentation order
%                   and with relevant initial positions, and frequncies filled in

daqreset
load(fullfile(direct, 'ExpStruct.mat'))

% create results dir if necessary
dirCont = dir(fullfile(direct,'results')); 

if length(dirCont) < 1 % if nothing is called results
    mkdir(direct, 'results')
else
    resDirInd = find([dirCont(:).isdir]);
    if ~resDirInd %not empty but not a dir
        mkdir(direc, 'results')
    end        
end


if nargin < 5
    random = 1;
end


% THIS SECTION SEEMS TO CAUSE COMM PROBLEMS
% %check that patterns and posFunc in expStruct match the SD card
% Panel_com('sync_sd_info') %updates the SD.mat file based on SD card in controller
% pause(2)
% load 'C:\Users\gruntmane\Documents\MATLAB\XmegaController_Matlab_V13\controller\SD.mat'
% sdPatNames = SD.pattern.pattNames;
% expStPatNames = {expStruct.pattern(:).name};
% sdPosFuncNames = SD.function.posFunctionName;
% expStPosFuncNames = {expStruct.posFunc(:).name};
% 
% patNamesCheck = prod(cellfun(@(x,y) strcmp(x,y(1:(end-4))), expStPatNames, sdPatNames));
% posFuncNameCheck = prod(cellfun(@(x,y) strcmp(x,y(1:(end-4))), expStPosFuncNames, sdPosFuncNames));
% 
% if ~patNamesCheck || ~posFuncNameCheck
%     error('SD card content does not agree with ExpStruct')
% end

numFreq = length(freqVec);


if random % makes sure that protocols are presented in a random order but within categories 
    allcat = [expStruct.protocol(protInds).Category];
    uniCat = unique(allcat);
    catCell = cell(1, length(uniCat));
    
    for ii=1:length(uniCat)
        tmpProts = protInds(allcat == uniCat(ii));
        tmpFreqs = arrayfun(@(x) ones(length(tmpProts), 1)*x, freqVec, 'uniformoutput', 0);
        tmpFreqs = vertcat(tmpFreqs{:});
        tmpPandF = vertcat(repmat(tmpProts, 1, numFreq), tmpFreqs');
        tmpShuff = randperm(length(tmpPandF));
        catCell{ii} = tmpPandF(:,tmpShuff);
    end
    allShuff = randperm(length(uniCat));
    randPandF = [catCell{allShuff}];
    expOrder = randPandF(1, :);
    expFreq = randPandF(2,:);
else
    expOrder = repmat(protInds, 1, numFreq);
    expFreq = arrayfun(@(x) ones(length(protInds), 1)*x, freqVec, 'uniformoutput', 0);
    expFreq = vertcat(expFreq{:})';
end

nProt = length(expOrder);

% setting up the session

ses = daq.createSession('ni');
ses.addAnalogInputChannel('Dev1',0:4, 'Voltage');
ses.Channels(1).InputType = 'SingleEnded';
ses.Channels(2).InputType = 'SingleEnded';
ses.Channels(3).InputType = 'SingleEnded';
ses.Channels(4).InputType = 'SingleEnded';
ses.Channels(5).InputType = 'SingleEnded';
ses.Rate = 10000;
disp(['session rate is: ', num2str(ses.Rate)])

expResults = struct;
expProts = expStruct.protocol(expOrder);

% setting up the waitbar to cancel function
wbh = waitbar(0,'1','Name','Presenting Protocol',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(wbh,'canceling',0)



for ii=1:nProt
    
    relInd = expOrder(ii);
    currProt = expStruct.protocol(relInd);
    
    if currProt.InitialPosition(1) == -1  % these values are define in the beginning of createProtocolScrip
        offset = [iniPos, 1];
    end
    
    expProts(ii).InitialPosition = offset;
    expProts(ii).FuncFreqX = expFreq(ii);
    expProts(ii).FuncFreqY = expFreq(ii);
    expProts(ii).ProtNum = relInd;
    
    expResults(ii).protocolNum = relInd;
    
    tempDur = calculateDuration(expStruct, relInd, expFreq(ii));
    expProts(ii).Duration = tempDur;
    
    ses.DurationInSeconds = tempDur; 
    
    Panel_com('set_pattern_id', currProt.PatternID);   
    pause(0.03)
    Panel_com('set_position', offset);  
    pause(0.03)
    Panel_com('set_mode', [4, 4]);
    pause(0.03)
    Panel_com('set_posfunc_id', currProt.PosFuncX); 
    pause(0.03)
    Panel_com('set_posfunc_id', currProt.PosFuncY); 
    pause(0.03)
    Panel_com('set_funcx_freq' ,  expFreq(ii));
    pause(0.03)
    Panel_com('set_funcy_freq' ,  expFreq(ii));
    pause(0.03)
    Panel_com('send_gain_bias', [0, 0, 0, 0]);
    pause(0.03)
     
    protTraits = expStruct.protocol(relInd).Traits;
    
    fprintf('Presenting protocol %d (%d of %d): %s\n', relInd, ii, nProt, protTraits) 
    
    waitbar(ii/nProt,wbh,sprintf('%d of %d: protocol %d',ii, nProt, relInd))
    
    Panel_com('start')
    data = ses.startForeground;
    expResults(ii).data = data;
    %in case experiment crashes
    save(fullfile(direct,'results',['tempDiffFreqResultsProtocol',num2str(relInd), '.mat']), 'data')
    Panel_com('stop')
    
    if getappdata(wbh,'canceling')
        break
    end
    
    
    
end

timeStamp = datestr(clock, 'yyyymmdd_HH-MM-SS');


save(fullfile(direct, 'results', ['diffFreqExpInfo', timeStamp,]), ...
     'expResults', 'expStruct', 'expProts')

delete(fullfile(direct, 'results','tempDiffFreqResults*.mat'))
Panel_com('all_off')

delete(wbh)

end
    


%%%%%%%%%%%%%%%%%%

function durInSec = calculateDuration(expStruct, pInd, relFreq)

% This subfunction calculates the new duration for the protocol based on
% its general index (pInd) and the new freq (relFreq)

postStimDur = 0.5; % same as in createProtocolsScript  

posFunInd = expStruct.protocol(pInd).PosFuncX(2);
posFunLength = expStruct.posFunc(posFunInd).length;
tempRemove = 1000 - relFreq * postStimDur;
durInSec = (posFunLength - tempRemove)/relFreq; % in secs



end







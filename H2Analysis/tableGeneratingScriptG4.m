

%% Moving bar 


fAppearInd = 1; % since posFunc starts from zero


relTab = protocolStruct.gratingTable;
span = relTab.span;
ort = relTab.orient;
wid = relTab.width;
frameFac = 1; % since in G4 it is handled with the position function

numStim = height(relTab);
barInds= zeros(numStim,3); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    
    relSp = span(ii);
    if rem(ort(ii),2) % to correct for diagonal movement
        relSp = 2*round(relSp/sqrt(2))+1;
    end

    fDis = relSp * frameFac + (wid(ii)-1) * frameFac - 1;
    cen = floor((fDis + fApp) /2); %since in fDis it is already gone
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
    barInds(ii,3) = cen + fAppearInd;
end
    
protocolStruct.gratingTable.firstFN = barInds(:,1);
protocolStruct.gratingTable.lastFN = barInds(:,2);
protocolStruct.gratingTable.cenFN = barInds(:,3);


clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen fAppearInd wid

protocolStruct.gratingTable

%% grating (diff pos) 


fAppearInd = 1; % since posFunc starts from zero


relTab = protocolStruct.gratingTable;
numFr = protocolStruct.inputParams.width * 2; %to generate full cycle
numCyc = protocolStruct.inputParams.numCyc;

numStim = height(relTab);
barInds= zeros(numStim,2); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    

    fDis = numFr * numCyc -1;
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
end
    
protocolStruct.gratingTable.firstFN = barInds(:,1);
protocolStruct.gratingTable.lastFN = barInds(:,2);


clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen fAppearInd wid
clear numCyc numFr

protocolStruct.gratingTable

%% Moving bar cont switch


fAppearInd = 1; % since posFunc starts from zero


relTab = protocolStruct.gratingTable;
span = relTab.span;
ort = relTab.orient;
wid = relTab.width;
swPos = relTab.switchPos;
frameFac = 1; % since in G4 it is handled with the position function

numStim = height(relTab);
barInds= zeros(numStim,4); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    
    relSp = span(ii);
    if rem(ort(ii),2) % to correct for diagonal movement
        relSp = 2*round(relSp/sqrt(2))+1;
    end

    fDis = relSp * frameFac + (wid(ii)-1) * frameFac - 1;
    cen = floor((fDis + fApp) /2); %since in fDis it is already gone
    swFrame = floor(relSp/2) + swPos(ii); % switchPos makes sure the switch corrects for width
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
    barInds(ii,3) = cen + fAppearInd;
    barInds(ii,4) = swFrame + fAppearInd;
end
    
protocolStruct.gratingTable.firstFN = barInds(:,1);
protocolStruct.gratingTable.lastFN = barInds(:,2);
protocolStruct.gratingTable.cenFN = barInds(:,3);
protocolStruct.gratingTable.switchFN = barInds(:,4);


clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen fAppearInd wid swPos swFrame

protocolStruct.gratingTable

%% Moving bar diff edge (limited orientation to 0:2:6 so not considered


% was frame 2 due to coding mistake (fixed for recordings after 21 Nov 18) 
fAppearInd = 1; % since posFunc starts from zero
% change back to 1

relTab = protocolStruct.gratingTable;
span = relTab.span;
wid = relTab.width;
frameFac = 1; % since in G4 it is handled with the position function

numStim = height(relTab);
barInds= zeros(numStim,3); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    relSp = span(ii);
    fDis = relSp * frameFac + (wid(ii)-1) * frameFac -1;
    cen = (fDis + fApp) /2; %since in fDis it is already gone
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
    barInds(ii,3) = cen + fAppearInd;
end
    
protocolStruct.gratingTable.firstFN = barInds(:,1);
protocolStruct.gratingTable.lastFN = barInds(:,2);
protocolStruct.gratingTable.cenFN = barInds(:,3);


clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen fAppearInd wid

protocolStruct.gratingTable

%% Moving edge diff pos 


fAppearInd = 1; % since posFunc starts from zero


relTab = protocolStruct.gratingTable;
span = relTab.span;
ort = relTab.orient;

frameFac = 1; % since in G4 it is handled with the position function

numStim = height(relTab);
barInds= zeros(numStim,3); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    
    relSp = span(ii);
    if rem(ort(ii),2) % to correct for diagonal movement
        relSp = 2*round(relSp/sqrt(2))+1;
    end

    fDis = relSp * frameFac -1;
    cen = floor((fDis + fApp) /2); %since in fDis it is already gone
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
    barInds(ii,3) = cen + fAppearInd;
end
    
protocolStruct.gratingTable.firstFN = barInds(:,1);
protocolStruct.gratingTable.lastFN = barInds(:,2);
protocolStruct.gratingTable.cenFN = barInds(:,3);


clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen fAppearInd wid

protocolStruct.gratingTable

%% cenSurr 

fAppearInd = 1; % since posFunc starts from zero

relTab = protocolStruct.gratingTable;
numStim = height(relTab);
barInds= zeros(numStim,2); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    
    fDis = fApp + 1;
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
end
    
protocolStruct.gratingTable.firstFN = barInds(:,1);
protocolStruct.gratingTable.lastFN = barInds(:,2);


clear barInds fDis fApp numStim  relTab ii fAppearInd 


protocolStruct.gratingTable

%% Flickerbar 


flickApp = 1;

relTab = protocolStruct.gratingTable;

totTime = protocolStruct.inputParams.flickerDur;
flickStV = protocolStruct.inputParams.flickerStart;

numStim = height(relTab);
onInds= cell(numStim,1); 
offInds = onInds;
disInd = zeros(numStim, 1);
stepsPerCyc = 2;

for ii=1:numStim
    
    tempNumCyc = floor(totTime/relTab.cycDur(ii));
    disInd(ii) = tempNumCyc * stepsPerCyc + flickApp -1; % since it is the last frame and not where it disappears
    if flickStV ==0
        offInds{ii} = flickApp:stepsPerCyc:disInd(ii)-1;
        onInds{ii} = flickApp+1:stepsPerCyc:disInd(ii);
    else
        offInds{ii} = flickApp+1:stepsPerCyc:disInd(ii);
        onInds{ii} = flickApp:stepsPerCyc:disInd(ii)-1;
    end
    
    
end



protocolStruct.gratingTable.onFN = onInds;
protocolStruct.gratingTable.offFN = offInds;
protocolStruct.gratingTable.firstFN = ones(numStim,1)*flickApp;
protocolStruct.gratingTable.lastFN = disInd;


clear relTab numStim ii tempNumStep *Inds frameFac totTime *App stepDurUnit disInd
clear stepsPerCyc flickFV tempNumCyc


protocolStruct.gratingTable


%% singleBar table


stimDur = protocolStruct.gratingTable.stimDur;
appear = ones(size(stimDur)) * 20;
%appear = ones(size(stimDur)) * 12;
disappear = zeros(size(stimDur));

disFrame = [20, 21, 23, 27]+1;
%disFrame = [12, 13, 15, 19]+1;

uValStimDur = unique(stimDur);

for ii=1:length(uValStimDur)
    
    disappear(stimDur == uValStimDur(ii)) = disFrame(ii);
    
end

protocolStruct.gratingTable.appear = appear;
protocolStruct.gratingTable.disappear = disappear;


clear disFrame uValStimDur stimDur appear disappear ii

%% Minimal motion 1

% not correct

fAppearInd = 1;


simFrame = protocolStruct.inputParams.stepDur / stepDurUnit;


relTab = protocolStruct.gratingTable;
span = abs(relTab.FBPos - relTab.SBPos);
frameFac = ceil(relTab.timeDiff); % in case there is a zero there (2 bars flashing together)

numStim = height(relTab);
barInds= zeros(numStim,4); % for appearance and disappearance of both bars

for ii=1:numStim
    
    if frameFac(ii) == 0
        fApp = 0;
        fDis = simFrame;
        sApp = 0;
        sDis = simFrame;
    else
        
        if span(ii) == 0
            fApp = frameFac(ii);
            sApp = fApp;
        else
            fApp = 0;    
            sApp = span(ii) * frameFac(ii);
        end
        
        fDis = fApp + frameFac(ii);
        
        sDis = sApp + frameFac(ii);
        
    end
    
    
   
    barInds(ii, :)  = [fApp, fDis, sApp, sDis] + fAppearInd;
end
    


protocolStruct.gratingTable.fAppear = barInds(:,1);
protocolStruct.gratingTable.fDisappear = barInds(:,2);
protocolStruct.gratingTable.sAppear = barInds(:,3);
protocolStruct.gratingTable.sDisappear = barInds(:,4);

clear barInds *Dis *App numStim spDiff relTab fAppearInd stepDurUnit simFrame ii frameFac



%% Minimal motion non speed corrected 



assert(all(protocolStruct.gratingTable.FBStat == 0), 'table was not corrected for FBstat=1')

fAppearInd = 1;

simFrame = protocolStruct.inputParams.stepDur;


relTab = protocolStruct.gratingTable;
span = abs(relTab.FBPos - relTab.SBPos);
fbV = relTab.FBVal;
sbV = relTab.SBVal;
frameFac = ceil(relTab.timeDiff); % in case there is zero (2 bars flashing together)

numStim = height(relTab);
barInds= zeros(numStim,4); % for appearance and disappearance of both bars

for ii=1:numStim
    
    if frameFac(ii) == 0
        fApp = 0;
        fDis = simFrame;
        sApp = 0;
        sDis = simFrame;
    else
        
        if span(ii) == 0 && fbV(ii) == sbV(ii) % since if they are different they both flash
            fApp = frameFac(ii);
            sApp = fApp;
        else
            fApp = 0;    
            sApp = frameFac(ii);
        end
        
        fDis = fApp + frameFac(ii);
        
        sDis = sApp + frameFac(ii);
        
    end
    
    
   
    barInds(ii, :)  = [fApp, fDis, sApp, sDis] + fAppearInd;
end


    
protocolStruct.gratingTable.fAppear = barInds(:,1);
protocolStruct.gratingTable.fDisappear = barInds(:,2);
protocolStruct.gratingTable.sAppear = barInds(:,3);
protocolStruct.gratingTable.sDisappear = barInds(:,4);

clear barInds *Dis *App numStim spDiff relTab fAppearInd stepDurUnit simFrame ii frameFac

protocolStruct.gratingTable



%%  minMot inhibition



fAppearInd = 20;
% fAppearInd = 12;
stepDurUnit = 0.02;
simFrame = protocolStruct.inputParams.stepDur / stepDurUnit;


relTab = protocolStruct.gratingTable;
span = abs(relTab.FBPos - relTab.SBPos);
frameFac = relTab.timeDiff/ stepDurUnit;
fbStat = relTab.FBStat;

numStim = height(relTab);
barInds= zeros(numStim,4); % for appearance and disappearance of both bars

for ii=1:numStim
    
    if frameFac(ii) == 0
        fApp = 0;
        fDis = simFrame;
        sApp = 0;
        sDis = simFrame;
    else
        
        if span(ii) == 0
            fApp = frameFac(ii);
            sApp = fApp;
        else
            fApp = 0;    
            sApp = span(ii) * frameFac(ii);
        end
        
        if fbStat(ii) == 0
            fDis = fApp + frameFac(ii);
        elseif fbStat(ii) == 1
            fDis = sApp + frameFac(ii);
        end
        
        sDis = sApp + frameFac(ii);
        
    end
    
    
   
    barInds(ii, :)  = [fApp, fDis, sApp, sDis] + fAppearInd;
end
    
protocolStruct.gratingTable.fAppear = barInds(:,1);
protocolStruct.gratingTable.fDisappear = barInds(:,2);
protocolStruct.gratingTable.sAppear = barInds(:,3);
protocolStruct.gratingTable.sDisappear = barInds(:,4);

clear barInds *Dis *App numStim spDiff relTab fAppearInd stepDurUnit simFrame ii frameFac fbStat


%% MinMoving Bar


fAppearInd = 12;
stepDurUnit = 0.02;



relTab = protocolStruct.gratingTable;
span = abs(relTab.startPos - relTab.stopPos);
frameFac = relTab.stepDur/ stepDurUnit;

numStim = height(relTab);
barInds= zeros(numStim,3); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    fDis = (span(ii)+1) * frameFac(ii);
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
    barInds(ii,3) = frameFac(ii);
end
    
protocolStruct.gratingTable.appear = barInds(:,1);
protocolStruct.gratingTable.disappear = barInds(:,2);
protocolStruct.gratingTable.framePerStep = barInds(:,3);

clear barInds *Dis *App numStim spDiff relTab ii frameFac fAppearInd stepDurUnit





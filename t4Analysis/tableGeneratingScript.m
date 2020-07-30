%% universal parameters

universalDur = [0.02, 0.04, 0.08, 0.16, 0.32];


%% singleBar table (and grating phase)


stimDur = protocolStruct.gratingTable.stimDur;
appear = ones(size(stimDur)) * 20;
%appear = ones(size(stimDur)) * 12;
disappear = zeros(size(stimDur));

disFrame = [20, 21, 23, 27, 35]+1;
%disFrame = [12, 13, 15, 19]+1;

uValStimDur = unique(stimDur);

for ii=1:length(uValStimDur)
    
    durInd = find(universalDur == uValStimDur(ii));
    
    disappear(stimDur == uValStimDur(ii)) = disFrame(durInd);
    
end

protocolStruct.gratingTable.appear = appear;
protocolStruct.gratingTable.disappear = disappear;


clear disFrame uValStimDur stimDur appear disappear ii uValStimDur disFrame disappear appear stimDur durInd


% uncomment for grating phase
% valSum  = protocolStruct.gratingTable.FBval + protocolStruct.gratingTable.SBVal;
% % switching indices so that 1+2=3 for the construction of the grating
% grtComb = round(valSum *2); 
% grtComb(grtComb == 2) = 4; 
% grtComb(grtComb == 3) = 2;
% grtComb(grtComb == 4) = 3;
% protocolStruct.gratingTable.grtComb = grtComb; 

protocolStruct.gratingTable

%% Minimal motion 1

tabString = protocolStruct.gratingTable.Properties.Description; 
widInd = strfind(tabString, 'Wid:') + 4; % since it gives index of first chr 
barWids = arrayfun(@(x) str2double(tabString(x)), widInd);
assert(diff(barWids) == 0, 'Deals only with same width bars')
relWid = barWids(1); 

fAppearInd = 20;
% fAppearInd = 12;
stepDurUnit = 0.02;
simFrame = protocolStruct.inputParams.stepDur / stepDurUnit;


relTab = protocolStruct.gratingTable;
span = abs(relTab.FBPos - relTab.SBPos);
frameFac = relTab.timeDiff/ stepDurUnit;

fbv = relTab.FBVal;
sbv = relTab.SBVal;
fbStat = relTab.FBStat;

numStim = height(relTab);
barInds= zeros(numStim,5); % for appearance and disappearance of both bars

for ii=1:numStim
    
    
    
    if frameFac(ii) == 0
        fApp = 0;
        fDis = simFrame;
        sApp = 0;
        sDis = simFrame;
        
        if fbv(ii) ~= sbv(ii) && span(ii) == 0% since in this case first bar will not be presented
            fApp = nan;
            fDis = nan;
        end
        
    else
        
        if span(ii) == 0
            if fbv(ii) == sbv(ii) % since if value is the same only the second is presented
                fApp = frameFac(ii);
                sApp = fApp;
            else % have to change for new version of minMot (recordings after Feb19)
                sApp = frameFac(ii);
                fApp = 0; 
            end
                
        else
            fApp = 0;    
            sApp = span(ii) * frameFac(ii);
        end
        
        fDis = fApp + frameFac(ii);
        sDis = sApp + frameFac(ii);
        
        if fbStat(ii) == 0
            fDis = fApp + frameFac(ii);
        elseif fbStat(ii) == 1
            fDis = sApp + frameFac(ii);
        end
        
    end
    
    if span(ii) == 0 
        overlap = -1;
    elseif span(ii) < relWid
        overlap = 1;
    else
        overlap = 0;
    end
    
   
    barInds(ii, 1:4)  = [fApp, fDis, sApp, sDis] + fAppearInd;
    barInds(ii, 5)  = overlap;
end
    
protocolStruct.gratingTable.fAppear = barInds(:,1);
protocolStruct.gratingTable.fDisappear = barInds(:,2);
protocolStruct.gratingTable.sAppear = barInds(:,3);
protocolStruct.gratingTable.sDisappear = barInds(:,4);
protocolStruct.gratingTable.overlap = barInds(:,5);

clear barInds *Dis *App numStim spDiff relTab fAppearInd stepDurUnit simFrame ii frameFac
clear sbv fbv overlap tabString relWid barWids widInd span fbStat

protocolStruct.gratingTable

%% Minimal motion non-speed corrected


fAppearInd = 20;
stepDurUnit = 0.02;

tabString = protocolStruct.gratingTable.Properties.Description; 
corrInd = strfind(tabString, 'Corr:') + 5; % since it gives index of first chr 
corrFlag = arrayfun(@(x) str2double(tabString(x)), corrInd);
assert(corrFlag == 0, 'Deals only with non speed corrected protocols')

widInd = strfind(tabString, 'Wid:') + 4; % to add the overlap column
barWids = arrayfun(@(x) str2double(tabString(x)), widInd);
assert(diff(barWids) == 0, 'Deals only with same width bars')
relWid = barWids(1); 


relTab = protocolStruct.gratingTable;
span = abs(relTab.FBPos - relTab.SBPos);
frameFac = relTab.timeDiff/ stepDurUnit;

fbv = relTab.FBVal;
sbv = relTab.SBVal;
fbStat = relTab.FBStat;

simFrame = protocolStruct.inputParams.stepDur / stepDurUnit;

numStim = height(relTab);
barInds= zeros(numStim,5); % for appearance and disappearance of both bars

for ii=1:numStim
    
    if frameFac(ii) == 0
        fApp = 0;
        fDis = simFrame;
        sApp = 0;
        sDis = simFrame;
        
        if fbv(ii) ~= sbv(ii) && span(ii) == 0% since in this case first bar will not be presented
            fApp = nan;
            fDis = nan;
        end
        
    else
        
        sApp = frameFac(ii);
        sDis = sApp + frameFac(ii);
        
        if span(ii) == 0
            if fbv(ii) == sbv(ii)
                fApp = sApp;
                fDis = sDis;
            else
                fApp = 0;
                fDis = frameFac(ii);
            end
                
        else
            fApp = 0;
            fDis = fApp + frameFac(ii);
        end
        
        if fbStat(ii) == 0
            fDis = fApp + frameFac(ii);
        elseif fbStat(ii) == 1
            fDis = sApp + frameFac(ii);
        end
        
    end
    
    if span(ii) == 0 
        overlap = -1;
    elseif span(ii) < relWid
        overlap = 1;
    else
        overlap = 0;
    end
    
   
    barInds(ii, 1:4)  = [fApp, fDis, sApp, sDis] + fAppearInd;
    barInds(ii, 5) = overlap; 
end
    
protocolStruct.gratingTable.fAppear = barInds(:,1);
protocolStruct.gratingTable.fDisappear = barInds(:,2);
protocolStruct.gratingTable.sAppear = barInds(:,3);
protocolStruct.gratingTable.sDisappear = barInds(:,4);
protocolStruct.gratingTable.overlap = barInds(:,5);

clear barInds *Dis *App numStim spDiff relTab fAppearInd stepDurUnit simFrame ii frameFac
clear overlap tabString relWid barWids span fbv sbv widInd corrInd corrFlag fbStat

protocolStruct.gratingTable

%% Minimal motion non speed corrected (dark) - old

% tabString = protocolStruct.gratingTable.Properties.Description; 
% widInd = strfind(tabString, 'Wid:') + 4; % since it gives index of first chr 
% barWids = arrayfun(@(x) str2double(tabString(x)), widInd);
% assert(diff(barWids) == 0, 'Deals only with same width bars')
% relWid = barWids(1); 
% 
% relTab = protocolStruct.gratingTable;
% span = abs(relTab.FBPos - relTab.SBPos);
% fbv = relTab.FBVal;
% sbv = relTab.SBVal;
% 
% 
% fAppearInd = 20;
% % fAppearInd = 12;
% stepDurUnit = 0.02;
% simFrame = protocolStruct.inputParams.stepDur / stepDurUnit;
% 
% 
% relTab = protocolStruct.gratingTable;
% span = abs(relTab.FBPos - relTab.SBPos);
% frameFac = relTab.timeDiff/ stepDurUnit;
% 
% numStim = height(relTab);
% barInds= zeros(numStim,5); % for appearance and disappearance of both bars
% 
% for ii=1:numStim
%     
%     if frameFac(ii) == 0
%         fApp = 0;
%         fDis = simFrame;
%         sApp = 0;
%         sDis = simFrame;
%     else
%         
%         sApp = 2 * frameFac(ii);
%         sDis = sApp + frameFac(ii);
%         
%         if span(ii) == 0
%             if fbv(ii) == sbv(ii)
%                 fApp = sApp;
%                 fDis = sDis;
%             else
%                 fApp = nan;
%                 fDis = nan;
%             end
%                 
%         else
%             fApp = 0;
%             fDis = fApp + frameFac(ii);
%         end
%     end
%     
%     if span(ii) == 0 
%         overlap = -1;
%     elseif span(ii) < relWid
%         overlap = 1;
%     else
%         overlap = 0;
%     end
%     
%    
%     barInds(ii, 1:4)  = [fApp, fDis, sApp, sDis] + fAppearInd;
%     barInds(ii, 5) = overlap; 
% end
%     
% protocolStruct.gratingTable.fAppear = barInds(:,1);
% protocolStruct.gratingTable.fDisappear = barInds(:,2);
% protocolStruct.gratingTable.sAppear = barInds(:,3);
% protocolStruct.gratingTable.sDisappear = barInds(:,4);
% protocolStruct.gratingTable.overlap = barInds(:,5);
% 
% clear barInds *Dis *App numStim spDiff relTab fAppearInd stepDurUnit simFrame ii frameFac
% clear overlap tabString relWid barWids span
% 
% protocolStruct.gratingTable


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

%% Flickerbar 


offApp = 20;
%offApp = 12;
stepDurUnit = 0.02;



relTab = protocolStruct.gratingTable;
frameFac = (relTab.cycDur/ stepDurUnit); %since it is full cycle (D & B) 
totTime = protocolStruct.inputParams.flickerDur;
flStVal = protocolStruct.inputParams.flickerStart;


numStim = height(relTab);
onInds= cell(numStim,1); 
offInds = onInds;
disInd = zeros(numStim, 1);


for ii=1:numStim
    
    tempNumStep = totTime/relTab.cycDur(ii) - 1;
    offInds{ii} = offApp:frameFac(ii):offApp+tempNumStep*frameFac(ii);
    onInds{ii} = offApp+frameFac(ii)/2:frameFac(ii):offApp+frameFac(ii)/2+tempNumStep*frameFac(ii);
    disInd(ii) = onInds{ii}(end)+frameFac(ii)/2;
    
end

% if starts bright need to flip indices
if flStVal
    tempInds = onInds; 
    onInds = offInds;
    offInds = tempInds; 
end

protocolStruct.gratingTable.onAppear = onInds;
protocolStruct.gratingTable.offAppear = offInds;
protocolStruct.gratingTable.appear = ones(numStim,1)*offApp;
protocolStruct.gratingTable.disappear = disInd;

clear relTab numStim ii tempNumStep *Inds frameFac totTime *App stepDurUnit disInd flStVal

protocolStruct.gratingTable

%% Moving bar 


fAppearInd = 20;
% fAppearInd = 12;
stepDurUnit = 0.02;



relTab = protocolStruct.gratingTable;
span = relTab.span;
ort = relTab.orient;
wid = relTab.width;
frameFac = relTab.stepDur/ stepDurUnit;

numStim = height(relTab);
barInds= zeros(numStim,5); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    
    relSp = span(ii);
    if rem(ort(ii),2) % to correct for diagonal movement
        relSp = 2*round(relSp/sqrt(2))+1;
    end
    
    cen = floor(relSp/2) * frameFac(ii); % since disappearing right before center is appearing in center
    fDis = relSp * frameFac(ii) + (wid(ii)-1) * frameFac(ii);
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
    barInds(ii,3) = frameFac(ii);
    barInds(ii,4) = cen + fAppearInd;
    barInds(ii,5) = barInds(ii,4) + wid(ii) * frameFac(ii);
end
    
protocolStruct.gratingTable.appear = barInds(:,1);
protocolStruct.gratingTable.disappear = barInds(:,2);
protocolStruct.gratingTable.framePerStep = barInds(:,3);
protocolStruct.gratingTable.leadCenter = barInds(:,4);
protocolStruct.gratingTable.trailCenter = barInds(:,5);


clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen 
clear fAppearInd wid span ort stepDurUnit 

protocolStruct.gratingTable

%% Moving Edge 


fAppearInd = 20;
% fAppearInd = 12;
stepDurUnit = 0.02;



relTab = protocolStruct.gratingTable;
span = relTab.span;
ort = relTab.orient;

frameFac = relTab.stepDur/ stepDurUnit;

numStim = height(relTab);
barInds= zeros(numStim,4); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    
    relSp = span(ii);
    if rem(ort(ii),2) % to correct for diagonal movement
        relSp = 2*round(relSp/sqrt(2))+1;
    end
    
    cen = floor(relSp/2) * frameFac(ii); % since disappearing right before center is appearing in center
    fDis = relSp * frameFac(ii);
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
    barInds(ii,3) = frameFac(ii);
    barInds(ii,4) = cen + fAppearInd;
end
    
protocolStruct.gratingTable.appear = barInds(:,1);
protocolStruct.gratingTable.disappear = barInds(:,2);
protocolStruct.gratingTable.framePerStep = barInds(:,3);
protocolStruct.gratingTable.center = barInds(:,4);


clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen fAppearInd wid ort stepDurUnit

valSum  = protocolStruct.gratingTable.edgeVal + protocolStruct.gratingTable.bkgdVal;
% switching indices so that 1+2=3 for the construction of the grating
grtComb = round(valSum *2); 
grtComb(grtComb == 2) = 4; 
grtComb(grtComb == 3) = 2;
grtComb(grtComb == 4) = 3;
protocolStruct.gratingTable.grtComb = grtComb; 


clear grtComb valSum grtStep

protocolStruct.gratingTable


%% moving grating

fAppearInd = 20;
stepDurUnit = 0.02;


numCyc = protocolStruct.inputParams.numCyc; 
grtSteps = 2*protocolStruct.inputParams.grtHWidth; % full cyc 
relTab = protocolStruct.gratingTable;



frameFac = relTab.stimDur/ stepDurUnit;

numStim = height(relTab);
barInds= zeros(numStim,3); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    
    fDis = numCyc * grtSteps * frameFac(ii);
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
    barInds(ii,3) = frameFac(ii);
end
    
protocolStruct.gratingTable.appear = barInds(:,1);
protocolStruct.gratingTable.disappear = barInds(:,2);
protocolStruct.gratingTable.framePerStep = barInds(:,3);


clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen fAppearInd wid numCyc


valSum  = protocolStruct.gratingTable.FBval + protocolStruct.gratingTable.SBVal;
% switching indices so that 1+2=3 for the construction of the grating
grtComb = round(valSum *2); 
grtComb(grtComb == 2) = 4; 
grtComb(grtComb == 3) = 2;
grtComb(grtComb == 4) = 3;
protocolStruct.gratingTable.grtComb = grtComb; 


clear grtComb valSum grtSteps

protocolStruct.gratingTable




%% Moving Figure Ground 

% works for wid 2 and 4 fpr fig and grd



fAppearInd = 20;
% fAppearInd = 12;
stepDurUnit = 0.02;



relTab = protocolStruct.gratingTable;
relParam = protocolStruct.inputParams;
posDif = relTab.posDiff; 
span = relParam.barSpan;
ort = relTab.orient;
fWid = relParam.figBarWid;
gWid = relParam.grdBarWid;

assert(fWid == 2 && gWid == 4, 'designed only for fWid =2 and gWid=4')

frameFac = relTab.stepDur/ stepDurUnit;

relVal = relTab.figureVal; 

numStim = height(relTab);
bkgInds= zeros(1,4); % for appearance and disappearance of both bars
figInds= zeros(1,4); 

barInds= zeros(numStim,9);


for ii=1:numStim
    
    figInds = zeros(1,4);
    bkgInds = figInds;
    
    if posDif(ii) == 0
       fgApp = frameFac(ii) * fWid; %since ground will appear only after fig
       effGWid = 2;
       posFac = 2;
    elseif posDif(ii) == 1
        effGWid = gWid;
        fgApp = 0;    
        posFac = 0;
    else
        effGWid = 2;
        fgApp = 0;  
        posFac = 0;
    end
    
    ffApp = frameFac(ii) * posDif(ii); 
    
    relSp = span;
    if rem(ort(ii),2) % to correct for diagonal movement
        relSp = 2*round(relSp/sqrt(2))+1;
    end
    
    cen = floor(relSp/2) * frameFac(ii); % since disappearing right before center is appearing in center
    
    fgDis = (relSp + posFac + effGWid-1) * frameFac(ii);
    ffDis = (relSp + fWid-1 + posDif(ii)) * frameFac(ii);
    
    bkgInds(1:2)  = [fgApp, fgDis] + fAppearInd;
    bkgInds(3) = cen + fgApp+ fAppearInd;
    bkgInds(4) = bkgInds(3) + effGWid * frameFac(ii);
    
    barInds(ii,9) = frameFac(ii);
    
    figInds(1:2)  = [ffApp, ffDis] + fAppearInd;
    figInds(3) = cen + ffApp + fAppearInd;
    figInds(4) = figInds(3) + fWid * frameFac(ii);
    
    if relVal(ii)
        barInds(ii, 1:4) = figInds; 
        barInds(ii, 5:8) = bkgInds; 
    else
        barInds(ii, 1:4) = bkgInds; 
        barInds(ii, 5:8) = figInds; 
    end
end
    

protocolStruct.gratingTable.onAppear = barInds(:,1);
protocolStruct.gratingTable.onDisappear = barInds(:,2);
protocolStruct.gratingTable.onLeadCenter = barInds(:,3);
protocolStruct.gratingTable.onTrailCenter = barInds(:,4);
protocolStruct.gratingTable.offAppear = barInds(:,5);
protocolStruct.gratingTable.offDisappear = barInds(:,6);
protocolStruct.gratingTable.offLeadCenter = barInds(:,7);
protocolStruct.gratingTable.offTrailCenter = barInds(:,8);
protocolStruct.gratingTable.framePerStep = barInds(:,9);

clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen fAppearInd *Wid
clear ffApp ffDis fgApp fgDis posDif posFac relVal relParam figInds bkgInds

protocolStruct.gratingTable

%% Moving bar shift


fAppearInd = 20;
% fAppearInd = 12;
stepDurUnit = 0.02;



relTab = protocolStruct.gratingTable;


poSeqLen = cellfun(@length, relTab.posSeq); 
wid = relTab.width;
frameFac = relTab.stepDur/ stepDurUnit;

numStim = height(relTab);
barInds= zeros(numStim,3); % for appearance and disappearance of both bars

for ii=1:numStim
    fApp = 0;    
    fDis = poSeqLen(ii) * frameFac(ii);
    
    barInds(ii, 1:2)  = [fApp, fDis] + fAppearInd;
    barInds(ii,3) = frameFac(ii);

end
    
protocolStruct.gratingTable.appear = barInds(:,1);
protocolStruct.gratingTable.disappear = barInds(:,2);
protocolStruct.gratingTable.framePerStep = barInds(:,3);


clear barInds fDis fApp numStim spDiff relTab ii frameFac relSp cen 
clear fAppearInd wid span ort stepDurUnit 

protocolStruct.gratingTable



%%

cn = 16; 
pNum = 1;

fileN = fullfile(t4NewCells(cn).dir, t4NewCells(cn).flicker(pNum).fileName);

load(fileN)

protocolStruct.gratingTable

%%

save(fileN, 'protocolStruct', '-v7.3')
disp('Saved!')

%%

clear protocolStruct* fileN






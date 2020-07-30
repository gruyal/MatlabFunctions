function expStruct = compileExpData(relDir)

% imports the TDMS files into matlab and combines them in a single
% structure
% relDir is the relevant directory for that experiment

tdmsFiles = dir(fullfile(relDir, '*.tdms'));
frameCh = floor(length(tdmsFiles)/2);

datCell = cell(frameCh, 2);

convTime = 10^-6;

for ii=1:length(tdmsFiles)
    
    tempFN = tdmsFiles(ii).name(1:end-5);
    chFlag = contains(tempFN, 'ADC');
    frFlag = contains(tempFN, 'Frame');
    
    if  chFlag || frFlag
        
        tempStruct = TDMS2struct(relDir,tempFN);
        gName = tempStruct.groupname{1}; 
        if chFlag
            relCh = str2double(gName(end))+1;
        else
            relCh = frameCh;
        end
        
        if strcmp(tempStruct.channame{1}, 'Time')
            relInd = 1;
        else
            relInd = 2;
        end
        
        tempDat = squeeze(tempStruct.chanvals(1,1,:));
        
        datCell{relCh, relInd} = tempDat;
        
    else
        
        tempStruct = TDMS_readTDMSFile(fullfile(relDir,[tempFN, '.tdms'])); % since TDMS2struct cant deal with this
        commTime = (tempStruct.data{end-2})';
        commType = (tempStruct.data{end-1})';
        commVal = (tempStruct.data{end})';
        firstTime = commTime(1);
        commTime = (commTime - firstTime)*convTime; %converts to sec
        commTable = table(commTime, commType, commVal);
    end
        
    
    
    
end

assert(unique(datCell{1,1} - datCell{2,1}) == 0, 'times are not equal across channels')

expStruct.commTable = commTable;

allDatC = datCell(1:end-1,2)';
allDatM = horzcat((datCell{1,1} - double(firstTime))* convTime, allDatC{:});

expStruct.chData = allDatM;

frameCell = datCell(end, :);
frameMat = [(frameCell{1} - double(firstTime)) * convTime, frameCell{2}];

expStruct.frameData = frameMat;




end


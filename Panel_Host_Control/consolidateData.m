function fullProtStruct = consolidateData(direct)

% function fullProtStruct = consolidateData(direct)
%
% This function loads the protocol structure in the given directory, takes
% the TDMS file names from it, and gets the relevant data from them. This
% is then saved as a single field in the loaded protocol structure. 
%
% INPUT
% direct - directory that contains both the TDMS files and a protocol structure
% with .stim.fileNames field, whose name starts with protocolStruct (the
% way runDumpProtocol saves it)
%
% OUTPUT
% fuuProtStruct - after the protocol structure is loaded it is saved with
% the relevant data (Analog inputs, and x Position).

% confirm that there is only one protocolStruct in the folder
pTest = dir(fullfile(direct, 'protocolStruct*.mat'));
assert(~isempty(pTest), 'no protocol structure in directory')
assert(length(pTest) == 1, 'More than one protocol structure in directory')


load(fullfile(direct, pTest.name));

for ii=1:length(protocolStruct.stim)
    
    fname = protocolStruct.stim(ii).fileName;
    if isempty(fname)
        fprintf('No file for stimulus %d \n', ii)
        continue
    end
    tdmsSt = TDMS_readTDMSFile(fullfile(direct, fname));
    inputChDat = extractADCDataFromTDMSStruct(tdmsSt);
    xPosDat = extractPositionFromTDMSStruct(tdmsSt);
    protocolStruct.stim(ii).data{1} = inputChDat;
    protocolStruct.stim(ii).data{2} = xPosDat;
    clear tdmsSt 
    fprintf('Exctraced data from stim %d of %d\n', ii, length(protocolStruct.stim))
    
end

save(fullfile(direct, pTest.name), 'protocolStruct')
fprintf('ProtocolStruct was saved with input data consolidated \n')
fullProtStruct = protocolStruct;    




end
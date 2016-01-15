function outStruct = generateDefaultPulseStruct

% function puStruct = generateDefaultPulseStruct
% 
% This function generates a default pulse structure to be used in
% makeAOVecForLED


puStruct.numPulse = {10, '', [1 inf]}; % number
puStruct.pulseWid = {10, '', [1 inf]}; % pulse width in ms
puStruct.ipi = {10, '', [0 inf]}; % inter pulse interval in ms
puStruct.firstPulse = {250, ''}; % Time of first pulase in ms
puStruct.padEndLen = {10000, '', [0 50000]}; % zeros to pad in the end of train


outStruct = StructDlg(puStruct, 'Pulse Struct');

end
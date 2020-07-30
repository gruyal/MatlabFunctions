
%% construct corrcetCorruptTable - Moving bar shift

% cell 3
cn = 3;
stim =   [50; 58; 70];
refRep = [ 2;  2;  2];
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MBS';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = table(cellNum, protType, stim, refRep); 

clear stim refRep protType cellNum cn

% cell 6
cn = 6;
stim = 54;
refRep = 2;
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MBS';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn

% cell 7
cn = 7;
stim =   [39; 40];
refRep = [ 3;  3];
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MBS';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn

% cell 9
cn = 9;
stim =   [22; 38];
refRep = [ 2;  1];
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MBS';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn


% cell 11
cn = 11;
stim =   [26; 33; 34; 37; 38; 42; 45; 46; 48];
refRep = [ 3;  1;  2;  2;  3;  1;  1;  1;  1];
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MBS';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn

% cell 12
cn = 12;
stim =   [27; 36; 43; 44; 46];
refRep = [ 1;  2;  2;  2;  2];
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MBS';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn


% cell 13
cn = 13;
stim =   40; 
refRep = 1; 
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MBS';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn


% cell 14
cn = 14;
stim =   [34; 38; 41]; 
refRep = [ 1;  1;  1]; 
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MBS';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn ii


%% Grating 

% cell 4
cn = 4;
stim =   [13; 15; 16; 23; 24]; 
refRep = [ 2;  1; -1; -3; -2]; 
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MG';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn ii


% cell 9
cn = 9;
stim =   [ 5; 7; 24]; 
refRep = [-2; 2; -2]; 
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MG';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn ii



% cell 11
cn = 11;
stim =   [ 7; 8; 14; 23]; 
refRep = [ 3; 2;  1;  2]; 
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MG';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn ii


% cell 12
cn = 12;
stim =   [ 10; 13; 14; 15; 24]; 
refRep = [ -3; -1; -2;  -2; 1]; % maybe 24 -3 (- is to remove that repeat  
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MG';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn ii


% cell 13
cn = 13;
stim =   13; 
refRep = 3;
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MG';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn ii


% cell 14
cn = 14;
stim =   [ 5; 7; 13; 14; 15; 23];
refRep = [-2; -2; -2; -2; 2; -3];
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MG';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn ii


% cell 16
cn = 16;
stim =   [ 3; 13; 23; 24];
refRep = [ 2; -3; -1; -3];
protType = cell(length(stim), 1); 
for ii=1:length(stim)
    protType{ii} = 'MG';
end
cellNum = ones(length(stim), 1) * cn;


corrCorruptTab = [corrCorruptTab; table(cellNum, protType, stim, refRep)]; 

clear stim refRep protType cellNum cn ii









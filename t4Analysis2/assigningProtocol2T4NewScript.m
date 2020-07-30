cd /Users/gruntmane/Documents/Research/ExpCodeandRes/panelController/PanelContExp/


%%

% dir that appear in both lists have 2 cells in them and one is useless

% uselessDir = { '20181221'; '20181231'; '20190102'; '20190103'; '20190107'; '20190116'; '20190123'};

% doubleCellDir = '20190123'; %second cell

relDirNames = { ...
                '20180914'; '20180915'; '20181203'; '20181204'; '20181205'; ...
                '20190102'; '20190103'; '20190118'; '20190121'; '20190122'; ...
                '20190122'; '20190125'; '20190129'; '20190131'; '20190201'; ...
                '20190207';
                };




%% 

cn=1;

    t4NewCells(cn).dir = relDirNames{cn};
    
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20180914_17-38/protocolStruct20180914_17-38.mat'; 
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20180914_17-45/protocolStruct20180914_17-45.mat'; 
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:1x5_Val:1'; % moved arena 
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20180914_17-47/protocolStruct20180914_17-47.mat'; 
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:3x3_Val:0,1';
    t4NewCells(cn).cenSurr(4).fileName = '/CenterSurroundProtocol20180914_17-49/protocolStruct20180914_17-49.mat'; 
    t4NewCells(cn).cenSurr(4).properties = 'Posmat:5x5_Val:1';
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20180914_17-53/protocolStruct20180914_17-53.mat';
    t4NewCells(cn).movBar(1).properties = 'W:2,4_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20180914_17-58/protocolStruct20180914_17-58.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20180914_18-01/protocolStruct20180914_18-01.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20180914_18-26/protocolStruct20180914_18-26.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0,1_trajNum:4'; 
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20180914_17-37/protocolStructAO20180914_17-37.mat';
    t4NewCells(cn).stepData(1).properties = 'old version';
    

%% 

cn=2;

    t4NewCells(cn).dir = relDirNames{cn};

    
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20180915_14-11/protocolStruct20180915_14-11.mat'; 
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20180915_14-26/protocolStruct20180915_14-26.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1'; % some more protocols in between - moved arena down
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20180915_14-28/protocolStruct20180915_14-28.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1';
    t4NewCells(cn).cenSurr(4).fileName = '/CenterSurroundProtocol20180915_14-31/protocolStruct20180915_14-31.mat';
    t4NewCells(cn).cenSurr(4).properties = 'Posmat:5x5_Val:1'; % moved a bit
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20180915_14-34/protocolStruct20180915_14-34.mat';
    t4NewCells(cn).movBar(1).properties = 'W:1,2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20180915_14-39/protocolStruct20180915_14-39.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20180915_14-42/protocolStruct20180915_14-42.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20180915_15-12/protocolStruct20180915_15-12.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    t4NewCells(cn).movBarShift(2).fileName = '/MovingBarDiagCorrShiftProtocol20180915_15-23/protocolStruct20180915_15-23.mat';
    t4NewCells(cn).movBarShift(2).properties = 'W:2,4_Dur:2,4_Val:0,1_offset:0_trajNum:12';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20180915_15-38/protocolStruct20180915_15-38.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:0,2,4_FB:B_SB:B';
    % not useful 
%     t4NewCells(cn).minMot(2).fileName = '/MinMotDiagCorrProtocol20180915_15-47/protocolStruct20180915_15-47.mat';
%     t4NewCells(cn).minMot(2).properties = 'W:4_Dur:2,4_FB:D_SB:B';
%     t4NewCells(cn).minMot(2).notes = 'was too noisy - aborted';
    
%%
    
cn=3;     
    
    t4NewCells(cn).dir = relDirNames{cn};

    
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20181203_17-09/protocolStruct20181203_17-09.mat'; 
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20181203_17-12/protocolStruct20181203_17-12.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1'; 
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20181203_17-15/protocolStruct20181203_17-15.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1';
    t4NewCells(cn).cenSurr(4).fileName = '/CenterSurroundProtocol20181203_17-18/protocolStruct20181203_17-18.mat';
    t4NewCells(cn).cenSurr(4).properties = 'Posmat:5x5_Val:1'; 
    t4NewCells(cn).cenSurr(4).notes= 'diff center, but previous better'; 
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20181203_17-21/protocolStruct20181203_17-21.mat';
    t4NewCells(cn).movBar(1).properties = 'W:2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20181203_17-23/protocolStruct20181203_17-23.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20181203_17-28/protocolStruct20181203_17-28.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20181203_17-59/protocolStruct20181203_17-59.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:1,2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20181203_18-11/protocolStruct20181203_18-11.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    
    t4NewCells(cn).movGrtWin(1).fileName = '/MovingGratingDiagCorrDiffWinProtocol20181203_18-20/protocolStruct20181203_18-20.mat';
    t4NewCells(cn).movGrtWin(1).properties = 'Dur:2,4_Win:9,15,21';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20181203_18-27/protocolStruct20181203_18-27.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    
%%

    
cn=4;     
    
    t4NewCells(cn).dir = relDirNames{cn};

    
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20181204_17-30/protocolStruct20181204_17-30.mat'; 
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20181204_17-35/protocolStruct20181204_17-35.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:1'; %checking center
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20181204_17-36/protocolStruct20181204_17-36.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:3x3_Val:0,1';
    t4NewCells(cn).cenSurr(4).fileName = '/CenterSurroundProtocol20181204_17-38/protocolStruct20181204_17-38.mat';
    t4NewCells(cn).cenSurr(4).properties = 'Posmat:5x5_Val:1'; 
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20181204_17-43/protocolStruct20181204_17-43.mat';
    t4NewCells(cn).movBar(1).properties = 'W:2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20181204_17-45/protocolStruct20181204_17-45.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20181204_17-49/protocolStruct20181204_17-49.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20181204_18-26/protocolStruct20181204_18-26.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20181204_18-20/protocolStruct20181204_18-20.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20181204_18-37/protocolStruct20181204_18-37.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20181204_17-28/protocolStructAO20181204_17-28.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20181204_17-30/protocolStructAO20181204_17-30.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:2,4,5_HypP:-1';
    t4NewCells(cn).stepData(3).fileName = '/stepData20181204_17-42/protocolStructAO20181204_17-42.mat';
    t4NewCells(cn).stepData(3).properties = 'Dur:2,4,5_HypP:-1';
    t4NewCells(cn).stepData(4).fileName = '/stepData20181204_18-17/protocolStructAO20181204_18-17.mat';
    t4NewCells(cn).stepData(4).properties = 'Dur:2,4,5_HypP:-1';
    t4NewCells(cn).stepData(5).fileName = '/stepData20181204_18-46/protocolStructAO20181204_18-46.mat';
    t4NewCells(cn).stepData(5).properties = 'Dur:2,4,5_HypP:-1';

    
    
    %%

    
cn=5;     
    
    t4NewCells(cn).dir = relDirNames{cn};
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20181205_16-30/protocolStructAO20181205_16-30.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20181205_16-33/protocolStructAO20181205_16-33.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:0,2,4,5_HypP:-1';
    
%%
    
cn=6;     
    
    t4NewCells(cn).dir = relDirNames{cn};
 
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190102_16-08/protocolStruct20190102_16-08.mat';
    t4NewCells(cn).cenSurr(1).properties = 'Posmat:3x3_Val:1'; %checking center
    t4NewCells(cn).cenSurr(1).notes = 'arena down';
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190102_16-09/protocolStruct20190102_16-09.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1';
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20190102_16-11/protocolStruct20190102_16-11.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1'; 
    t4NewCells(cn).cenSurr(3).notes = 'nicely centered';
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190102_16-14/protocolStruct20190102_16-14.mat';
    t4NewCells(cn).movBar(1).properties = 'W:2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190102_16-16/protocolStruct20190102_16-16.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190102_16-20/protocolStruct20190102_16-20.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    t4NewCells(cn).singleBar(3).fileName = '/SingleBarDiagCorrDiffWProtocol20190102_17-30/protocolStruct20190102_17-30.mat';
    t4NewCells(cn).singleBar(3).properties = 'W:2,4_Dur:2,4_Val:1_Posover:1';
    t4NewCells(cn).singleBar(3).notes = 'room lights on, used flashlight a bit';
    t4NewCells(cn).singleBar(4).fileName = '/SingleBarDiagCorrDiffWProtocol20190102_17-46/protocolStruct20190102_17-46.mat';
    t4NewCells(cn).singleBar(4).properties = 'W:2,4_Dur:2,4_Val:0,1_Posover:0';
    t4NewCells(cn).singleBar(4).notes = 'big flashlight pointing at fly';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190102_16-49/protocolStruct20190102_16-49.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0,1_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20190102_17-05/protocolStruct20190102_17-05.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    
    t4NewCells(cn).movGrtWin(1).fileName = '/MovingGratingDiagCorrDiffWinProtocol20190102_17-12/protocolStruct20190102_17-12.mat';
    t4NewCells(cn).movGrtWin(1).properties = 'Dur:2,4_Win:9,15,21';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20190102_17-17/protocolStruct20190102_17-17.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    t4NewCells(cn).movGrtPhase(1).notes= 'arena FU on contra side';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190102_15-49/protocolStructAO20190102_15-49.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20190102_17-25/protocolStructAO20190102_17-25.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:0,2,4,5_HypP:-1';
    t4NewCells(cn).stepData(3).fileName = '/stepData20190102_17-56/protocolStructAO20190102_17-56.mat';
    t4NewCells(cn).stepData(3).properties = 'Dur:0,2,4,5_HypP:-1';
    
    
    
%%
    
    
cn=7;

    t4NewCells(cn).dir = relDirNames{cn};
    
    % has more cenSurr since I had to move arena around
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190103_17-33/protocolStruct20190103_17-33.mat';
    t4NewCells(cn).cenSurr(1).properties = 'Posmat:3x3_Val:0,1';
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190103_17-35/protocolStruct20190103_17-35.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:5x5_Val:1';
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190103_17-39/protocolStruct20190103_17-39.mat';
    t4NewCells(cn).movBar(1).properties = 'W:1,2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190103_17-42/protocolStruct20190103_17-42.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190103_17-45/protocolStruct20190103_17-45.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    t4NewCells(cn).singleBar(2).notes = 'had to increase inj in the middle - resp filtered'; 
    t4NewCells(cn).singleBar(3).fileName = '/SingleBarDiagCorrDiffWProtocol20190103_18-20/protocolStruct20190103_18-20.mat';
    t4NewCells(cn).singleBar(3).properties = 'W:2_Dur:4_Val:1_Posover:1';
    t4NewCells(cn).singleBar(3).notes = 'repeated to check quality, was reduced - aborted'; 
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190103_18-12/protocolStruct20190103_18-12.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0,1_trajNum:4'; 

    
%%


cn=8;     
    
    t4NewCells(cn).dir = relDirNames{cn};
 
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190118_13-23/protocolStruct20190118_13-23.mat';
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190118_13-27/protocolStruct20190118_13-27.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1';
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20190118_13-29/protocolStruct20190118_13-29.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1'; 
    t4NewCells(cn).cenSurr(4).fileName = '/CenterSurroundProtocol20190118_13-32/protocolStruct20190118_13-32.mat';
    t4NewCells(cn).cenSurr(4).properties = 'Posmat:5x5_Val:1'; %rechecking center
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190118_13-35/protocolStruct20190118_13-35.mat';
    t4NewCells(cn).movBar(1).properties = 'W:1,2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190118_13-38/protocolStruct20190118_13-38.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190118_13-41/protocolStruct20190118_13-41.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190118_14-15/protocolStruct20190118_14-15.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20190118_14-09/protocolStruct20190118_14-09.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20190118_14-26/protocolStruct20190118_14-26.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190118_13-21/protocolStructAO20190118_13-21.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20190118_14-35/protocolStructAO20190118_14-35.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:0,2,4,5_HypP:-1';
    t4NewCells(cn).stepData(3).fileName = '/stepData20190118_14-39/protocolStructAO20190118_14-39.mat';
    t4NewCells(cn).stepData(3).properties = 'Dur:0,2,4,5_HypP:-1';
    t4NewCells(cn).stepData(3).notes = 'room lights on';
    
    
    
%%


cn=9;     
    
    t4NewCells(cn).dir = relDirNames{cn};
 
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190121_15-03/protocolStruct20190121_15-03.mat';
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190121_15-08/protocolStruct20190121_15-08.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1';
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20190121_15-11/protocolStruct20190121_15-11.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1'; 
    t4NewCells(cn).cenSurr(4).fileName = '/CenterSurroundProtocol20190121_15-13/protocolStruct20190121_15-13.mat';
    t4NewCells(cn).cenSurr(4).properties = 'Posmat:5x5_Val:1'; %different center
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190121_15-16/protocolStruct20190121_15-16.mat';
    t4NewCells(cn).movBar(1).properties = 'W:2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190121_15-19/protocolStruct20190121_15-19.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190121_15-22/protocolStruct20190121_15-22.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    t4NewCells(cn).singleBar(3).fileName = '/SingleBarDiagCorrDiffWProtocol20190121_16-28/protocolStruct20190121_16-28.mat';
    t4NewCells(cn).singleBar(3).properties = 'W:2,4_Dur:2,4_Val:0,1_Posover:1';
    t4NewCells(cn).singleBar(3).notes = 'room lights on';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190121_15-56/protocolStruct20190121_15-56.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20190121_15-50/protocolStruct20190121_15-50.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20190121_16-09/protocolStruct20190121_16-09.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    
    t4NewCells(cn).movGrtWin(1).fileName = '/MovingGratingDiagCorrDiffWinProtocol20190121_16-18/protocolStruct20190121_16-18.mat';
    t4NewCells(cn).movGrtWin(1).properties = 'Dur:2,4_Win:9,15,21';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190121_15-01/protocolStructAO20190121_15-01.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20190121_16-06/protocolStructAO20190121_16-06.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:0,2,4,5_HypP:-1';
    

 
 %%


cn=10;     
    
    t4NewCells(cn).dir = relDirNames{cn};
 
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190122_14-36/protocolStruct20190122_14-36.mat';
    t4NewCells(cn).cenSurr(1).notes  = 'after this moved arena all the way down'; 
    %more Center srround for non useful conditions
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190122_14-42/protocolStruct20190122_14-42.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1';
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20190122_14-45/protocolStruct20190122_14-45.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1'; 
    t4NewCells(cn).cenSurr(4).fileName = '/CenterSurroundProtocol20190122_14-48/protocolStruct20190122_14-48.mat';
    t4NewCells(cn).cenSurr(4).properties = 'Posmat:5x5_Val:1'; %different center
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190122_14-50/protocolStruct20190122_14-50.mat';
    t4NewCells(cn).movBar(1).properties = 'W:2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190122_14-55/protocolStruct20190122_14-55.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190122_14-57/protocolStruct20190122_14-57.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    t4NewCells(cn).singleBar(2).notes = 'noisy - aborted cell after';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190122_14-30/protocolStructAO20190122_14-30.mat';

    
%%


cn=11;     
    
    t4NewCells(cn).dir = relDirNames{cn};
 
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190122_16-37/protocolStruct20190122_16-37.mat';
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190122_16-41/protocolStruct20190122_16-41.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1';
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20190122_16-43/protocolStruct20190122_16-43.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1';
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190122_16-46/protocolStruct20190122_16-46.mat';
    t4NewCells(cn).movBar(1).properties = 'W:1,2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190122_16-49/protocolStruct20190122_16-49.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190122_16-52/protocolStruct20190122_16-52.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190122_17-21/protocolStruct20190122_17-21.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20190122_17-30/protocolStruct20190122_17-30.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20190122_17-37/protocolStruct20190122_17-37.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    
    t4NewCells(cn).movGrtWin(1).fileName = '/MovingGratingDiagCorrDiffWinProtocol20190122_17-46/protocolStruct20190122_17-46.mat';
    t4NewCells(cn).movGrtWin(1).properties = 'Dur:2,4_Win:9,15,21';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190122_16-35/protocolStructAO20190122_16-35.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20190122_17-35/protocolStructAO20190122_17-35.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:0,2,4,5_HypP:-1';
    
    
    
%%    

cn=12;     
    
% check RF is not cut by bottom edge of arena 

    t4NewCells(cn).dir = relDirNames{cn};
 
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190125_14-45/protocolStruct20190125_14-45.mat';
    t4NewCells(cn).cenSurr(1).notes = 'arena down'; 
    % arena malfunction in between - several cenSurr not compeltely sure 
    
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190125_14-58/protocolStruct20190125_14-58.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:5x5_Val:1';
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190125_15-00/protocolStruct20190125_15-00.mat';
    t4NewCells(cn).movBar(1).properties = 'W:1,2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190125_15-04/protocolStruct20190125_15-04.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    % there is another singlebar - but arena was cutting bar in half - not
    % useful 
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190125_15-26/protocolStruct20190125_15-26.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    % there is a 15:59 singlebar that was aborted due to arena malfunction
    t4NewCells(cn).singleBar(3).fileName = '/SingleBarDiagCorrDiffWProtocol20190125_16-05/protocolStruct20190125_16-05.mat';
    t4NewCells(cn).singleBar(3).properties = 'W:2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190125_15-52/protocolStruct20190125_15-52.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20190125_15-49/protocolStruct20190125_15-49.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    t4NewCells(cn).minMot(2).fileName = '/MinMotDiagCorrProtocol20190125_16-27/protocolStruct20190125_16-27.mat';
    t4NewCells(cn).minMot(2).properties = 'W:4_Dur:2,4_FB:D_SB:B - more positions';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20190125_16-15/protocolStruct20190125_16-15.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    
    t4NewCells(cn).movGrtWin(1).fileName = '/MovingGratingDiagCorrDiffWinProtocol20190125_16-22/protocolStruct20190125_16-22.mat';
    t4NewCells(cn).movGrtWin(1).properties = 'Dur:2,4_Win:9,15,21';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190125_14-39/protocolStructAO20190125_14-39.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20190125_15-19/protocolStructAO20190125_15-19.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:0,2,4,5_HypP:-2';
    t4NewCells(cn).stepData(3).fileName = '/stepData20190125_16-31/protocolStructAO20190125_16-31.mat';
    t4NewCells(cn).stepData(3).properties = 'Dur:0,2,4,5_HypP:-3';
    
    
    
    %%    

cn=13;

    t4NewCells(cn).dir = relDirNames{cn};
 
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190129_16-18/protocolStruct20190129_16-18.mat';
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190129_16-22/protocolStruct20190129_16-22.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1'; 
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20190129_16-24/protocolStruct20190129_16-24.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1';
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190129_16-27/protocolStruct20190129_16-27.mat';
    t4NewCells(cn).movBar(1).properties = 'W:1,2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190129_16-30/protocolStruct20190129_16-30.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190129_16-33/protocolStruct20190129_16-33.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190129_17-06/protocolStruct20190129_17-06.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20190129_17-00/protocolStruct20190129_17-00.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    t4NewCells(cn).minMot(2).fileName = '/MinMotDiagCorrProtocol20190129_17-32/protocolStruct20190129_17-32.mat';
    t4NewCells(cn).minMot(2).properties = 'W:4_Dur:2,4_FB:D,B_SB:B_FBStat:0,1';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20190129_17-17/protocolStruct20190129_17-17.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    
    t4NewCells(cn).movGrtWin(1).fileName = '/MovingGratingDiagCorrDiffWinProtocol20190129_17-27/protocolStruct20190129_17-27.mat';
    t4NewCells(cn).movGrtWin(1).properties = 'Dur:2,4_Win:9,15,21';
    
    t4NewCells(cn).flicker(1).fileName = '/FlickerBarDiagCorrProtocol20190129_17-38/protocolStruct20190129_17-38.mat';
    t4NewCells(cn).flicker(1).properties = 'Dur:2,4_W:4';
    t4NewCells(cn).flicker(1).notes = 'nice fast events';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190129_16-16/protocolStructAO20190129_16-16.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20190129_17-25/protocolStructAO20190129_17-25.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:0,2,4,5_HypP:-2';
   
  
    
    
    
    %%    

cn=14;

    t4NewCells(cn).dir = relDirNames{cn};
 
    
    % check that RF is not cut off at the top of arena
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190131_17-27/protocolStruct20190131_17-27.mat';
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190131_17-31/protocolStruct20190131_17-31.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1'; 
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20190131_17-33/protocolStruct20190131_17-33.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1';
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190131_17-37/protocolStruct20190131_17-37.mat';
    t4NewCells(cn).movBar(1).properties = 'W:1,2_Dur:3_Ori:8_Val:1';
    t4NewCells(cn).movBar(1).notes = 'different center';
    
    % there is another single bar with a diff center
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190131_17-43/protocolStruct20190131_17-43.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190131_17-46/protocolStruct20190131_17-46.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    t4NewCells(cn).singleBar(2).notes = 'check the edge resp - could be double cell again?!';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190131_18-22/protocolStruct20190131_18-22.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20190131_18-13/protocolStruct20190131_18-13.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    t4NewCells(cn).minMot(2).fileName = '/MinMotDiagCorrProtocol20190131_18-42/protocolStruct20190131_18-42.mat';
    t4NewCells(cn).minMot(2).properties = 'W:4_Dur:2,4_FB:D,B_SB:B_FBStat:0,1';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20190131_18-33/protocolStruct20190131_18-33.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    
    t4NewCells(cn).flicker(1).fileName = '/FlickerBarDiagCorrProtocol20190131_18-50/protocolStruct20190131_18-50.mat';
    t4NewCells(cn).flicker(1).properties = 'Dur:2,4_W:4';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190131_17-26/protocolStructAO20190131_17-26.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20190131_18-47/protocolStructAO20190131_18-47.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:0,2,4,5_HypP:-2';
   
    
    
     %%    

cn=15;

    t4NewCells(cn).dir = relDirNames{cn};
 
    
    % check that RF is not cut off at the top of arena
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190201_16-08/protocolStruct20190201_16-08.mat';
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190201_16-12/protocolStruct20190201_16-12.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:1x5_Val:1'; % arena down 
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20190201_16-13/protocolStruct20190201_16-13.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:3x3_Val:0,1'; 
    t4NewCells(cn).cenSurr(4).fileName = '/CenterSurroundProtocol20190201_16-15/protocolStruct20190201_16-15.mat';
    t4NewCells(cn).cenSurr(4).properties = 'Posmat:5x5_Val:1';
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190201_16-18/protocolStruct20190201_16-18.mat';
    t4NewCells(cn).movBar(1).properties = 'W:1,2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190201_16-21/protocolStruct20190201_16-21.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190201_16-23/protocolStruct20190201_16-23.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190201_16-57/protocolStruct20190201_16-57.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20190201_16-50/protocolStruct20190201_16-50.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190201_16-06/protocolStructAO20190201_16-06.mat';
   
    
    
    
     %%    

cn=16;

    t4NewCells(cn).dir = relDirNames{cn};
 
    t4NewCells(cn).cenSurr(1).fileName = '/CenterSurroundProtocol20190207_18-01/protocolStruct20190207_18-01.mat';
    % several more cenSurr in adjusting arena height - arena down 
    t4NewCells(cn).cenSurr(2).fileName = '/CenterSurroundProtocol20190207_18-09/protocolStruct20190207_18-09.mat';
    t4NewCells(cn).cenSurr(2).properties = 'Posmat:3x3_Val:0,1'; 
    t4NewCells(cn).cenSurr(3).fileName = '/CenterSurroundProtocol20190207_18-11/protocolStruct20190207_18-11.mat';
    t4NewCells(cn).cenSurr(3).properties = 'Posmat:5x5_Val:1';
    
    t4NewCells(cn).movBar(1).fileName = '/MovingBarDiagCorrProtocol20190207_18-13/protocolStruct20190207_18-13.mat';
    t4NewCells(cn).movBar(1).properties = 'W:1,2_Dur:3_Ori:8_Val:1';
    
    t4NewCells(cn).singleBar(1).fileName = '/SingleBarDiagCorrDiffWProtocol20190207_18-17/protocolStruct20190207_18-17.mat';
    t4NewCells(cn).singleBar(1).properties = 'W:2_Dur:4_Val:1_Posover:1'; %presinglebar
    t4NewCells(cn).singleBar(2).fileName = '/SingleBarDiagCorrDiffWProtocol20190207_18-20/protocolStruct20190207_18-20.mat';
    t4NewCells(cn).singleBar(2).properties = 'W:1,2,4_Dur:2,4_Val:0,1_Posover:1';
    
    t4NewCells(cn).movBarShift(1).fileName = '/MovingBarDiagCorrShiftProtocol20190207_18-52/protocolStruct20190207_18-52.mat';
    t4NewCells(cn).movBarShift(1).properties = 'W:1,2,4_Dur:2,4_Val:0,1_offset:0_trajNum:4';
    
    t4NewCells(cn).minMot(1).fileName = '/MinMotDiagCorrProtocol20190207_18-47/protocolStruct20190207_18-47.mat';
    t4NewCells(cn).minMot(1).properties = 'W:4_Dur:2,4_FB:D_SB:B';
    t4NewCells(cn).minMot(2).fileName = '/MinMotDiagCorrProtocol20190207_19-00/protocolStruct20190207_19-00.mat';
    t4NewCells(cn).minMot(2).properties = 'W:4_Dur:2,4_FB:D,B_SB:B_FBStat:0,1';
    
    t4NewCells(cn).movGrtPhase(1).fileName = '/MovingGratingDiagCorrDiffPhaseProtocol20190207_19-07/protocolStruct20190207_19-07.mat';
    t4NewCells(cn).movGrtPhase(1).properties = 'Dur:2,4_3Cont';
    t4NewCells(cn).movGrtPhase(1).notes = 'nice resp';
    
    t4NewCells(cn).movGrtWin(1).fileName = '/MovingGratingDiagCorrDiffWinProtocol20190207_19-22/protocolStruct20190207_19-22.mat';
    t4NewCells(cn).movGrtWin(1).properties = 'Dur:2,4_Win:9,15,21';
    t4NewCells(cn).movGrtWin(1).notes = 'nice resp';
    
    % another aborted flicker since start B
    t4NewCells(cn).flicker(1).fileName = '/FlickerBarDiagCorrProtocol20190207_19-17/protocolStruct20190207_19-17.mat';
    t4NewCells(cn).flicker(1).properties = 'Dur:2,4_W:2,4';
    
    t4NewCells(cn).stepData(1).fileName = '/stepData20190207_18-00/protocolStructAO20190207_18-00.mat';
    t4NewCells(cn).stepData(2).fileName = '/stepData20190207_19-21/protocolStructAO20190207_19-21.mat';
    t4NewCells(cn).stepData(2).properties = 'Dur:0,2,4,5_HypP:-1.5';
    
 %%
  
clear cn relDirNames

%%

allProtFields = fieldnames(t4NewCells); 

protByCell = zeros(length(t4NewCells), length(allProtFields)-1); 

for nn=2:length(allProtFields)
    
    for cc=1:length(t4NewCells)
        
        protByCell(cc, nn-1) = length(getfield(t4NewCells(cc), allProtFields{nn}));
        
    end
    
end
        

protByCellTab = array2table(protByCell, 'variablenames', allProtFields(2:end));

cellNo = (1:length(t4NewCells))';

protByCellTab = addvars(protByCellTab, cellNo, 'before', allProtFields{2});

clear nn cc cellNo protByCell allProtFields    


    
    





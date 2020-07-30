% Script version of Function_Maker_G4 with current GUI parameters
% (script saved in C:\matlabroot\G4\Scripts\)
%
% Save this script with a new filename to keep it from being overwritten

%% user-defined function parameters
pfnparam.type = 'pfn'; %number of frames in pattern
pfnparam.frames = 16; %number of frames in pattern
pfnparam.gs_val = 4; %brightness bits in pattern
pfnparam.section = { 'static' 'triangle' 'triangle' 'triangle' 'static' }; %static, sawtooth, traingle, sine, cosine, or square
pfnparam.dur = [ 0.25 2.5 2.5 2.5 0.25 ]; %section duration (in s)
pfnparam.val = [ 1 nan nan nan 1 ]; %function value for static sections
pfnparam.high = [ nan 16 16 16 nan ]; %high end of function range {for non-static sections}
pfnparam.low = [ nan 1 1 1 nan ]; %low end of function range {for non-static sections}
pfnparam.freq = [ 1 1 1 1 1 ]; %frequency of section {for non-static sections}
pfnparam.flip = [ 0 0 1 0 0 ]; %flip the range of values of function {for non-static sections}


%% generate function
func = Function_Maker_G4(pfnparam);


%% save function
save_dir = 'C:\matlabroot\G4\Position Functions\';
cd(save_dir);
funcfiles = ls('Function*.mat');
if isempty(funcfiles)
    pfnparam.ID = 1;
else
    pfnparam.ID = size(funcfiles,1)+1;
end
filename = ['Function_' num2str(pfnparam.ID, '%04d') '_G4.mat'];
save_function_G4(func, pfnparam, save_dir, filename);


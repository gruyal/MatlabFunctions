% Script to create directional tuning map
exp_name = 'exp_Dec17';

%% user-defined pattern parameters
% all angles/distances/sizes are in units of radians rather than degrees
% some parameters are only needed in certain cirumstances {specified by curly brace}
param.pattern_type = 'square grating'; %square grating, sine grating, edge, starfield, or Off/On
param.motion_type = 'rotation'; %rotation, translation, or expansion-contraction
param.pattern_fov = 'full-field'; %full-field or local
param.arena_pitch = deg2rad(0); %angle of arena pitch (0 = straight ahead, positive values = pitched up)
param.gs_val = 4; %bits of intensity value (1 or 4)
param.levels = [0 8 8]; %brightness level of [1st bar (in grating) or advancing edge, 2nd bar or receding edge, background (mask)]
param.pole_coord = [0 -pi/2]; %location of pattern pole [longitude, lattitude] {for pattern_fov=full-field}
param.motion_angle = deg2rad(0); %angle of rotation (0=rightward motion, positive values rotate the direction clockwise) {fov=local}
param.spat_freq = deg2rad(360); %spatial angle (in radians) before pattern repeats {for gratings and edge}
param.step_size = deg2rad(1.875); %amount of motion per frame (in radians) {for type~=off/on}
param.duty_cycle = 4.1667; %percent of spat_freq taken up by first bar {for square gratings} 
param.num_dots = 500; %number of dots in star-field {for type=starfield}
param.dot_radius = 0.02182; %radius of dots (in radians) {for starfield}
param.dot_size = 'static'; %static or distance-relative {for starfield}
param.dot_occ = 'closest'; %how occluding dots are drawn (closest, sum, or mean) {for starfield}
param.dot_level = 0; %0 = dot brightness set to 1st level; 1 and 2 = random brightness (0-1st; 0 or 1st) {for starfield}
param.snap_dots = 0; %1 if apparent dot locations should be rounded to the nearest pixel {for starfield}
param.sa_mask = deg2rad([0 0 180 0]); %location, size, and direction of solid angle mask [longitude, lattitude, solid_angle, out/in]
param.long_lat_mask = [-pi pi -pi/2 pi/2 0]; %coordinates of lattitude/longitude mask [min-long, max-long, min-lat, max-lattitude, out/in]
param.aa_samples = 15; %# of samples taken to calculate the brightness of each pixel (1 or 15 suggested)
param.aa_poles = 1; %1=anti-aliases the poles of rotation/translation grating/edge stimuli by matching them to the duty cycle
param.back_frame = 0; %1=adds a frame (frame 1) uniformly at background (mask) level
param.flip_right = 0; %1=left-right flips the right half of the pattern
param.phase_shift = 0; %shifts the starting frame of pattern


%% user-defined position function parameters

% rotFac = 4;

pfnparam.type = 'pfn'; %number of frames in pattern
pfnparam.frames = 17; %number of frames in pattern
pfnparam.gs_val = 4; %brightness bits in pattern
pfnparam.section = {'static', 'sawtooth', 'sawtooth', 'sawtooth', 'sawtooth', 'static'}; %static, sawtooth, triangle, sine, cosine, or square
pfnparam.dur = [0.5, 2.5, 2.5, 2.5, 2.5, 0.5]; %section duration (in s)
pfnparam.val = [ 1 nan nan nan nan 1 ]; %function value for static sections
pfnparam.high = [ nan 17 17 17 17 nan ]; %high end of function range {for non-static sections}
pfnparam.low = [ nan 2 2 2 2 nan ]; %low end of function range {for non-static sections}
pfnparam.freq = [ nan 9 9 9 9 nan ]; %frequency of section {for non-static sections}
pfnparam.flip = [ nan 0 1 0 1 nan ]; %flip the range of values of function {for non-static sections}



%% user-defined AO function parameters
% afnparam.type = 'afn'; %number of frames in pattern
% afnparam.section = { 'static' 'static' }; %static, sawtooth, traingle, sine, cosine, or square
% afnparam.dur = [ 0.1 4.9 ]; %section duration (in s)
% afnparam.val = [ 10 1 ]; %function value for static sections
% afnparam.high = [ 8 10 ]; %high end of function range {for non-static sections}
% afnparam.low = [ 1 2 ]; %low end of function range {for non-static sections}
% afnparam.freq = [ 1 9 ]; %frequency of section {for non-static sections}
% afnparam.flip = [ 0 0 ]; %flip the range of values of function {for non-static sections}

param.ID = 0;

%create experiment folders
exp_folder = create_exp_dir_G4(exp_name);

    
%% make stripe
param.ID = param.ID + 1;

%make and save pattern
[Pats, ~, ~] = Motion_Maker_G4(param);
param.stretch = zeros(size(Pats,3),1);
save_pattern_G4(Pats, param, [exp_folder '\Patterns'], ['Pattern_' num2str(param.ID, '%04d') '_G4.mat']);
% make and save position function
func = Function_Maker_G4(pfnparam); pfnparam.ID = param.ID;
save_function_G4(func, pfnparam, [exp_folder '\Functions'], ['Function_' num2str(param.ID, '%04d') '_G4.mat']);
% make and save analog output function
% func = Function_Maker_G4(afnparam); afnparam.ID = param.ID;
% save_function_G4(func, afnparam, [exp_folder '\Analog Output Functions'], ['FunctionAO_' num2str(param.ID, '%04d') '_G4.mat']);
    

%% make other grating stimuli
param.levels = [15 0 8];
param.spat_freq = deg2rad(30);
param.duty_cycle = 50;
param.back_frame = 1;

for speed = [1 2 4 8]
    for levels = [1 3 15]
        param.ID = param.ID + 1;
        pfnparam.freq = [ nan speed speed speed speed nan ];
        param.levels = [levels 0 8];

        %make and save pattern
        [Pats, ~, ~] = Motion_Maker_G4(param);
        param.stretch = zeros(size(Pats,3),1);
        save_pattern_G4(Pats, param, [exp_folder '\Patterns'], ['Pattern_' num2str(param.ID, '%04d') '_G4.mat']);
        % make and save position function
        func = Function_Maker_G4(pfnparam); pfnparam.ID = param.ID;
        crFunc = correctPositionFunction(func, pfnparam);
        
        if levels == 1 %checking correction
            plot(func)
            hold on 
            plot(crFunc)
            pause
            close(gcf)
        end
        
        save_function_G4(crFunc, pfnparam, [exp_folder '\Functions'], ['Function_' num2str(param.ID, '%04d') '_G4.mat']);
        % make and save analog output function
        % func = Function_Maker_G4(afnparam); afnparam.ID = param.ID;
        % save_function_G4(func, afnparam, [exp_folder '\Analog Output Functions'], ['FunctionAO_' num2str(param.ID, '%04d') '_G4.mat']);
    end
end


%% finalize experiment folder
create_currentExp(exp_folder)
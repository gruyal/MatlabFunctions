%% User-defined experiment conditions
experiment_name = 'exp_Dec17'; %  flightdemo_yawlift   pitchslip   rollthrust   'Speed_Tuning_2.5and5degSS_30degSF'
num_reps = 3;
randomize = 1; %1=yes, 0=no
exp_mode = 1; %1=position function, 2=constant rate, 3=position change, 4=CL, 5=CL+bias, 6=CL+OL
inter_type = 1; %0=no intertrial, 1=inter-trial closed loop with pattern 1
fly_name = 'DM_12';
inter_trial_duration = 8;
CL_gain = -20;
CL_offset = 0;


pre_start = 1;
stripe_ID = 1; % pattern ID for stripe fixation 
opto_ID = 7; % pattern ID for optomotor (during adjustement period)
fixDur = 20;
optoDur = 10; 

%% Load configuration
experiment_folder = ['C:\matlabroot\G4\Experiments\' experiment_name];
load([experiment_folder '\currentExp.mat']);
trial_duration = currentExp.trialDuration; %seconds
num_conditions = currentExp.pattern.num_patterns;
if ~exist(fullfile(experiment_folder,'Log Files'),'dir')
    mkdir(experiment_folder,'Log Files');
end

%% check if log files already present
if exist([experiment_folder '\Log Files\*'],'file')
    fprintf('unsorted log files present in save folder, remove before restarting experiment\n');
    return
end
if exist([experiment_folder '\Results\' fly_name],'dir')
    fprintf('Results folder already exists with that fly name\n');
    return
end

%% create .mat file of experiment order
skip_first = 0;
if inter_type==1
    skip_first = 1;
end
num_conditions = num_conditions - skip_first;

%set experiment order
if randomize == 1
    exp_order = NaN(num_reps,num_conditions);
    for rep_ind = 1:num_reps
        exp_order(rep_ind,:) = randperm(num_conditions);
    end
else
    exp_order = repmat(1:num_conditions,num_reps,1);
end
exp_order = exp_order + skip_first;


connectHost;
Panel_com('change_root_directory', experiment_folder);
exp_seconds = num_reps*num_conditions*(trial_duration+inter_trial_duration);
fprintf(['Estimated experiment duration: ' num2str(exp_seconds/60) ' minutes\n']);
save([experiment_folder '\Log Files\exp_order.mat'],'exp_order')
    

%% start experiment


if pre_start==1
    
    f = figure('Visible','off', 'position', [1065, 794, 159, 48]);
    f.MenuBar = 'none'; 
    f.ToolBar = 'none'; 

    btn = uicontrol('Style', 'pushbutton', 'String', 'Start Expariment',...
            'Position', [20 20 100 20],...
            'Callback', @pbCallback); 

    f.Visible = 'on';
    
    while ishandle(f)
        Panel_com('set_control_mode', 4);
        Panel_com('set_gain_bias', [CL_gain CL_offset]);
        Panel_com('set_pattern_id', stripe_ID);
        Panel_com('set_position_x',185)
        Panel_com('start_display', fixDur*10); %duration expected in 100ms units
        pause(fixDur)
        Panel_com('stop_display');
        
        if ishandle(f) % just to add another breakpoint

            Panel_com('set_control_mode', 1);
            Panel_com('set_pattern_id', opto_ID);
            Panel_com('set_pattern_func_id', opto_ID);
            Panel_com('start_display', optoDur * 10);
            pause(optoDur)
            Panel_com('stop_display');
            
        end
        
%         start = input('enter "Y" to start experiment (or "N" to continue stripe fixation): ','s');
%         
%         if strncmpi(start,'y',1)==0
%             break
%         end
    end
end





%%
% Panel_com('set_active_ao_channels', dec2bin(1,4)); %4 bits represent the active channels (e.g. 0001 = ch0, 0110 = ch2&3)
Panel_com('start_log');
pause(inter_trial_duration+1);


wbH = waitbar(0,'Repeat:  Stimulus:','Name','Progress',...
              'CreateCancelBtn',...
              'setappdata(gcbf,''canceling'',1)');
setappdata(wbH,'canceling',0)

totNum = num_reps * num_conditions;
count = 0;
tic
for r = 1:num_reps
    for c = 1:num_conditions
        cond_ID = exp_order(r,c);
        count = count+1;
        
        waitbar(count/totNum,wbH,sprintf('Repeat: %d Stimulus: %d',r,c))
        
        %trial portion
        Panel_com('set_control_mode', exp_mode);
        Panel_com('set_pattern_id', cond_ID);
        Panel_com('set_pattern_func_id', cond_ID);
%         Panel_com('set_ao_function_id',[1, cond_ID]); %argument(1) specifies channel: 0-3
%         fprintf(['Rep ' num2str(r) ' of ' num2str(num_reps) ', cond ' num2str(c) ' of ' num2str(num_conditions) ...
%             ': ' strjoin(currentExp.pattern.pattNames(cond_ID)) '\n']);
        Panel_com('start_display', trial_duration*10); %duration expected in 100ms units
        pause(trial_duration)
        Panel_com('stop_display');
        
        
        %intertrial portion
        if inter_type == 1
            Panel_com('set_control_mode', 4);
            Panel_com('set_gain_bias', [CL_gain CL_offset]);
            Panel_com('set_pattern_id', stripe_ID);
            Panel_com('set_position_x',185);
            %Panel_com('set_ao_function_id',[1, cond_ID]);
            Panel_com('start_display', inter_trial_duration*10);
            pause(inter_trial_duration);
            Panel_com('stop_display');
        end
        
        if getappdata(wbH,'canceling')
            break
        end
        
    end
    
    if getappdata(wbH,'canceling')
        break
    end
    
end

Panel_com('stop_display');
%rename/move results folder
pause(0.5);
Panel_com('stop_log');
pause(1);
movefile([experiment_folder '\Log Files\*'],fullfile(experiment_folder,'Results',fly_name));
disconnectHost;
delete(wbH)
toc
disp('finished');


%%

function pbCallback(btn, eventdata, handles)
 % hObject    handle to pushbutton1 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 disp('Starting Experiment');
 close(gcf);
 end
function alignFly

% function alignFly
% this function runs a set experiment with fixation and optomotor stimuli
% that are fixed into the function. It does not log the responses, only
% used for proper alignment

relDir = 'C:\matlabroot\G4_Display_Tools\G4_Protocol_Designer'; % experimental directory
fixPat = 3; %fixation pattern in this experiment
optPat = 1; % optomotor pattern
optFun = 2; % optomotor position function

% add check if it is already connected
connectHost;

Panel_com('change_root_directory', relDir);
inptStr = 'n';

while strcmpi('y', inptStr)

    % fixation 
    Panel_com('set_control_mode', 4);
    Panel_com('set_pattern_id', fixPat);

    Panel_com('start_display', 2000); %duration expected in 100ms units
    pause % wait for feedback


    % optomotor

    Panel_com('set_control_mode', 1);
    Panel_com('set_pattern_id', optPat);
    Panel_com('set_pattern_func_id', optFun);

    Panel_com('start_display', 2000);
    pause
    
    
    % repeats if requested
    prompt = 'Do you want to repeat the alignment? Y/N : ';
    inptStr = input(prompt,'s');
    
end


end
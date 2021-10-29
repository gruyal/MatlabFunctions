clear

% gratings stimuli
frame_dur = 80;
gr_spatial_period = 8;
gr_spatial_off = 4;
gr_spatial_bgr = gr_spatial_period-gr_spatial_off;
num_periods = 2;
gr_num_cycles = 4;
gr_num_pos = gr_spatial_period*num_periods;

single_period = [-1*ones(1,gr_spatial_off), 0*ones(1,gr_spatial_bgr)];
% single_period1 = repmat([-1*ones(1,gr_spatial_off), 0*ones(1,gr_spatial_off)], 1, gr_spatial_period/4);
single_frame = repmat(single_period,1,num_periods);
one_cycle =...
    toeplitz([single_frame(1) fliplr(single_frame(2:end))], single_frame);

single_frame1 = [single_frame(1) fliplr(single_frame(2:end))];
one_cycle1 = nan(length(single_frame1));
for ii=1:length(single_frame1)
    one_cycle1(ii,:) = circshift(single_frame1, (ii-1)*gr_spatial_off);
end


stim = repmat(one_cycle,gr_num_cycles,1);
num_time_steps = size(stim,1);
gr_pos = -floor(gr_num_pos/2):floor(gr_num_pos/2);
gr_pos = gr_pos(1:gr_num_pos);

rev_stim = stim;
rev_stim(1:2:num_time_steps,:) = -stim(1:2:num_time_steps,:);

time = (0:(num_time_steps-1))*frame_dur;

figure(1)
clf
subplot(2,2,1)
imagesc(gr_pos,time,stim)
caxis([-1,1])
title('standard PD')
subplot(2,2,2)
imagesc(gr_pos,time,fliplr(stim))
caxis([-1,1])
title('standard ND')
subplot(2,2,3)
imagesc(gr_pos,time,rev_stim)
caxis([-1,1])
title('reverse phi PD')
subplot(2,2,4)
imagesc(gr_pos,time,fliplr(rev_stim))
caxis([-1,1])
title('reverse phi ND')

%--------------------------------------------------------------
%run dynamics
frameR = 1;
modTimeVec = (0:frameR:(time(end)+frame_dur+30))'; %30: delay

load /Users/gruntmane/'Dropbox (HHMI)'/revPhiPaper/forSandro/relParamsT5.mat

% relParamsT5(21:28) = 0;

[V_pd, fe_pd, fe2_pd, fi_pd, fi2_pd,...
    stimOnset_pd, stepD_pd, modVal_pd, posVec_pd] =...
    run_stim(stim, time, gr_pos, modTimeVec, relParamsT5);

[V_nd, fe_nd, fe2_nd, fi_nd, fi2_nd,...
    stimOnset_nd, stepD_nd, modVal_nd, posVec_nd] =...
    run_stim(fliplr(stim), time, gr_pos, modTimeVec, relParamsT5);

[rev_V_pd, rev_fe_pd, rev_fe2_pd, rev_fi_pd, rev_fi2_pd,...
    rev_stimOnset_pd, rev_stepD_pd, rev_modVal_pd, rev_posVec_pd] =...
    run_stim(rev_stim, time, gr_pos, modTimeVec, relParamsT5);

[rev_V_nd, rev_fe_nd, rev_fe2_nd, rev_fi_nd, rev_fi2_nd,...
    rev_stimOnset_nd, rev_stepD_nd, rev_modVal_nd, rev_posVec_nd] =...
    run_stim(fliplr(rev_stim), time, gr_pos, modTimeVec, relParamsT5);

figure(3)
clf
subplot(2,1,1)
plot(modTimeVec, V_pd, 'linewidth', 2, 'color', 'k')
hold on
plot(modTimeVec, V_nd, 'linewidth', 2, 'color', 'r')
% plot(modTimeVec, fe2_nd, 'linewidth', 2, 'color', 'r', 'linestyle', '--')
% plot(modTimeVec, fi2_nd, 'linewidth', 2, 'color', 'r', 'linestyle', '--')
% plot(modTimeVec, fe2_pd, 'linewidth', 2, 'color', 'k', 'linestyle', '--')
% plot(modTimeVec, fi2_pd, 'linewidth', 2, 'color', 'k', 'linestyle', '--')
legend({'pd','nd'})
title('standard')
set(gca, 'ylim', [-10, 20])
subplot(2,1,2)
plot(modTimeVec, rev_V_pd, 'linewidth', 2, 'color', 'k')
hold on
plot(modTimeVec, rev_V_nd, 'linewidth', 2, 'color', 'r')
title('reverse phi')
set(gca, 'ylim', [-5, 10])


%--------------------------------------------------------------
%generate inputs from stimuli and call t5_simple_wrap_v2
function [V, fe, fe2, fi, fi2, stimOnset, stepD, modVal, posVec] =...
    run_stim(stim, time, gr_pos, modTimeVec, relParamsT5)
frame_dur = time(2)-time(1);
frameR = modTimeVec(2) - modTimeVec(1);
stimOnset = [];
posVec = {};
stepD = {};
modVal = {};
k = 1;
for i = 1:length(gr_pos)
    stim_single_pos = stim(:,i);
    last_stim_val = 0;
    dur_stim = 0;
    in_stim = 0;
    for j = 1:length(stim_single_pos)
        if(j==1) % first time step
            if(stim_single_pos(j)~=0) %if not background start stim
                stimOnset(k) =  time(j);
                posVec{k} = gr_pos(i);
                if(stim_single_pos(j)==1)
                    modVal{k} = 1;
                else
                    modVal{k} = 0;
                end
                last_stim_val = stim_single_pos(j);
                last_stim_time = time(j);
                dur_stim = frame_dur;
                in_stim = 1;
            end
        end
        if(j>1 && j~=length(stim_single_pos)) %intermediate time
            if(stim_single_pos(j)==last_stim_val) %same luminance as before
                if(in_stim) %and during stim, then increase stim dur
                    dur_stim = dur_stim + frame_dur;
                end
            end
            if(stim_single_pos(j)~=last_stim_val) %change of luminance
                if(in_stim) % and during stim then end stim
                    stepD{k} = dur_stim;
                    dur_stim = 0;
                    in_stim = 0;
                    last_stim_val = 0;
                    k = k + 1; %new index for next stim
                    if(stim_single_pos(j)~=0) %if ~ backgr, start new stim
                        stimOnset(k) =  time(j);
                        posVec{k} = gr_pos(i);
                        if(stim_single_pos(j)==1)
                            modVal{k} = 1;
                        else
                            modVal{k} = 0;
                        end
                        last_stim_val = stim_single_pos(j);
                        last_stim_time = time(j);
                        dur_stim = frame_dur;
                        in_stim = 1;
                    end
                else % (not in stim) then start stim
                    stimOnset(k) =  time(j);
                    posVec{k} = gr_pos(i);
                    if(stim_single_pos(j)==1)
                        modVal{k} = 1;
                    else
                        modVal{k} = 0;
                    end
                    last_stim_val = stim_single_pos(j);
                    last_stim_time = time(j);
                    dur_stim = frame_dur;
                    in_stim = 1;
                end
            end
        end
        if(j==length(stim_single_pos)) %last time
            if(in_stim && stim_single_pos(j)==last_stim_val)
                % and during stim then end stim
                stepD{k} = dur_stim + frame_dur;
                k = k + 1;
            else
                if(stim_single_pos(j)~=last_stim_val) %change of luminance
                    if(in_stim) % and during stim then end stim
                        stepD{k} = dur_stim;
                        dur_stim = 0;
                        in_stim = 0;
                        last_stim_val = 0;
                        k = k + 1; %new index for next stim
                        if(stim_single_pos(j)~=0)
                            %if ~ backgr, start and end new stim
                            stimOnset(k) =  time(j);
                            posVec{k} = gr_pos(i);
                            if(stim_single_pos(j)==1)
                                modVal{k} = 1;
                            else
                                modVal{k} = 0;
                            end
                            stepD{k} = frame_dur;
                            k = k + 1;
                        end
                    else % (not in stim) then start and end stim
                        stimOnset(k) =  time(j);
                        posVec{k} = gr_pos(i);
                        if(stim_single_pos(j)==1)
                            modVal{k} = 1;
                        else
                            modVal{k} = 0;
                        end
                        stepD{k} = frame_dur;
                        k = k + 1;
                    end   
                end
            end
        end
    end
end

stimInds = cell(1,length(stimOnset));
for i=1:length(stimOnset)
    tmp = zeros(size(modTimeVec));
    tmp(stimOnset(i)+1)=1;
    stimInds{i} = tmp;
end

[V, fe, fe2, fi, fi2] =...
    t5_simple_wrap_v2(relParamsT5,modVal,...
    stepD, posVec,modTimeVec,frameR,stimInds);

end




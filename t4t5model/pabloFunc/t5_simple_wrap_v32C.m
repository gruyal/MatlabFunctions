function [V, fe, fe2, fi, fi2] = t5_simple_wrap_v32C(...
    params,val_in,dur_in,pos_vect_in,time,fr,stimIdx_in,...
    jumpValOnset,jumpValOffset,ff_on,ff_off)
%{
T5 Simple Model

Output
    V:
      t x 1 vector with a simulation of T5 response to stimulus for a time
            vector t

Input
    params:
     a 1 x 29 vector of parameters in order [mu, sig, A, Tr, Td, me,
     be, bi, Ti]. Each parameter is a block fo values for [E, I, E2, I2].
                ie) [muE, muI, muE2, muI2, sigE, sigI, ..., Ti]
    val:
       Luminance of stimulus. 0 for dark relative to baseline, 1 for
                bright. paired cell with all stim_inputs
    dur:
       duration each step, in ms. Usually 40 or 160.
       paired cell with all stim_inputs
    pos_vect:
        a vector indicating the Receptive Field-centered positioins
        to stimulate, in integers, indexed from leading edge.
         paired cell with all stim_inputs
    time:
      tx1 vector with entry values of time course, in ms. If the
                to simulate 1000ms sampled with a framerate of 0.5 ms,
                [0:0.5:1000].
    fr:         framerate, in ms.
    stimIdx:
        logical of size t x 1, indexing when the next flash should
        be. sum(stimIdx) == length(pos_vect).
        paired cell with all stim_inputs
    %}
    stim_curr          = 0;
    % initialize a counter to finish loops when all stims are done
    
    stim_tot           = length(val_in);
    % finish repeating when we've done all stims
    
    ge_mat             = zeros(stim_tot,length(time));
    fe_mat             = ge_mat;
    fe2_mat            = ge_mat;
    gi_mat             = zeros(stim_tot,length(time));
    fi_mat             = gi_mat;
    fi2_mat            = gi_mat;
    % preallocate to store conductances to each stimulus
    
    params_to_keep = who;
    params_to_keep{end+1} = 'params_to_keep';
    while stim_curr < stim_tot
        clearvars('-except', params_to_keep{:})
        stim_curr = stim_curr+1; %go through each stim cond
        
        spfr_data.val      = val_in{stim_curr};
        %spfr_data.width    = width_in{stim_curr};
        spfr_data.stimDur  = dur_in{stim_curr};
        spfr_data.pos_vect = pos_vect_in{stim_curr};
        spfr_data.time     = time;
        spfr_data.fr       = fr;
        spfr_data.stimIdx  = stimIdx_in{stim_curr};
        spfr_data.jumpOn = jumpValOnset{stim_curr};
        spfr_data.jumpOff = jumpValOffset{stim_curr};
        if(jumpValOnset{stim_curr}==2)
            spfr_data.jumpOn = ff_on;
        end
        if(jumpValOffset{stim_curr}==2)
            spfr_data.jumpOff = ff_off;
        end
        
        if spfr_data.val == 0
            %all params stored in "rows" of 4,
            mu_p  = [params(8:9), 0, 0];
            sig_p = [params(10:11), 1, 1];
            amp_p = [params(6:7), 0, 0];
            tr_p  = params([4,5,4,5])*10;
            td_p  = params([2,3,2,3])*10;
            me_p  = 0.01;
            mi_p  = 0.01;
            be_p  = 0;
            bi_p  = 0;
            ti_p  = 1;
        else
            mu_p  = params([3:4,1:2]);
            %if NC, flip order of all param sets
            %(seconds always fed into delta)
            
            sig_p = [1,1,1,1];
            amp_p = [0,0,0,0];
            tr_p  = params([4,5,4,5])*10;
            td_p  = params([2,3,2,3])*10;
            me_p  = 0.01;
            mi_p  = 0.01;
            be_p  = 0;
            bi_p  = 0;
            ti_p  = 1;
        end
        
        p.Vr = 0;
        p.Ve = 65;
        p.Vi = -9;
        
        %stim parameters
        %width = spfr_data.width-1;
        %p.width = spfr_data.width;
        p.stimDur = spfr_data.stimDur;
        pos_vect = spfr_data.pos_vect;
        d = 30/spfr_data.fr; 
        %set manual delay to replicate delay in
        %transmission from photoreceptors to T4/5 response
        
        tmp = find(spfr_data.stimIdx);
        stimIdx = false(size(spfr_data.stimIdx));
        stimIdx(tmp+d) = true;
        stim_time = fix(spfr_data.time(stimIdx));
        
        M = 1:(length(spfr_data.pos_vect));
        x = min(pos_vect):max(pos_vect);
        %x =  [x0(1)+(-width:-1), x0]; 
        
        Ie = zeros(length(stim_time),length(M));
        Ii = zeros(length(stim_time),length(M));
        
        %pos_vect defines leadingedge of stimulated position.
        %x defines valid
        %positions in window
        for i = 1:length(pos_vect)
            %idx = x <= pos_vect(i) & x >= pos_vect(i)-width;
            idx = x <= pos_vect(i) & x >= pos_vect(i);
            Ie(i,idx) = 1;
            Ii(i,idx) = 1;
        end
        
        p.Ie = Ie;
        p.Ii = Ii;
        
        %excitatory (classic)
        p.Tre =  tr_p(1);
        p.Tde =  td_p(1);
        Ae    = amp_p(1);
        mue   =  mu_p(1);
        sige  = sig_p(1);
        
        %inhibitory (classic)
        p.Tri =  tr_p(2);
        p.Tdi =  td_p(2);
        Ai    = amp_p(2);
        mui   =  mu_p(2);
        sigi  = sig_p(2);
        
        %excitatory2 (delta, from PC)
        p.Tre2 =  tr_p(3);
        p.Tde2 =  td_p(3);
        Ae2    = amp_p(3);
        mue2   =  mu_p(3);
        sige2  = sig_p(3);
        
        %inhibitory2 (delta, from PC)
        p.Tri2 =  tr_p(4);
        p.Tdi2 =  td_p(4);
        Ai2    = amp_p(4);
        mui2   =  mu_p(4);
        sigi2  = sig_p(4);
        
        %time
        effDur = spfr_data.stimDur;
        
        %both contrasts are linear
        p.be2 = max(0,effDur*me_p  + be_p);
        p.bi2 = max(0,effDur*mi_p + bi_p);
        
        p.Tid = ti_p;
        
        % Model
        p.ae  =  Ae.*exp(-(x-mue).^2  / (2*sige^2) );
        p.ai  =  Ai.*exp(-(x-mui).^2  / (2*sigi^2) );
        p.ae2 = Ae2.*exp(-(x-mue2).^2 / (2*sige2^2) );
        p.ai2 = Ai2.*exp(-(x-mui2).^2 / (2*sigi2^2) );
        
        p.stim_time = stim_time;
        t_ind = find([stimIdx;1]);
        [~,s_ind] =...
            min(abs(spfr_data.time - (spfr_data.time(stimIdx) + effDur)'));
        s_ind = s_ind';
        
        p.stim_time = stim_time;
        
        %populate a matrix with this dynamic for each stim
        fe = zeros(size(p.Ie,2),size(spfr_data.time,1));
        fi = zeros(size(p.Ii,2),size(spfr_data.time,1));
        fe2 = zeros(size(p.Ii,2),size(spfr_data.time,1));
        fi2 = zeros(size(p.Ii,2),size(spfr_data.time,1));
        
        t_tot = spfr_data.time;
        t_on  = t_ind(1);
        t_off = s_ind(1);
        
        fe_tmp  = fc_sol(t_tot, t_on, t_off, p.Tre, p.Tde);
        fi_tmp  = fc_sol(t_tot, t_on, t_off, p.Tri, p.Tdi);
        
        fe2_tmp =...
            fd_sol(t_tot, t_on, t_off, p.Tid, p.Tre2, p.Tde2, p.be2);
        fi2_tmp =...
            fd_sol(t_tot, t_on, t_off, p.Tid, p.Tri2, p.Tdi2, p.bi2);
        
        %after solving for conductances at each time,
        %we can now shift back
        %to non-delayed frame to fill vector
        %{
         %S: this is unused (stimIdx has been shifted)
         t_ind = t_ind-d;
         s_ind = s_ind-d;
        %}
        
        for i = 1:size(Ie,1)
            fe(logical(p.Ie(i,:)),:)  = repmat(fe_tmp,sum(p.Ie(i,:)),1);
            fi(logical(p.Ii(i,:)),:)  = repmat(fi_tmp,sum(p.Ii(i,:)),1);
            fe2(logical(p.Ii(i,:)),:)  = repmat(fe2_tmp,sum(p.Ii(i,:)),1);
            fi2(logical(p.Ii(i,:)),:)  = repmat(fi2_tmp,sum(p.Ii(i,:)),1);
        end
        
        fe  =  fe(:,1:size(spfr_data.time,1));
        fi  =  fi(:,1:size(spfr_data.time,1));
        fe2 = fe2(:,1:size(spfr_data.time,1));
        fi2 = fi2(:,1:size(spfr_data.time,1));
        
        %caluclate conductances and steady state voltage
        
        ge_mat(stim_curr,:)  = p.ae*fe*spfr_data.jumpOn+...
            p.ae2*fe2*spfr_data.jumpOff;
        %ge_mat(stim_curr,:)  = p.ae*fe*1+...
        %    p.ae2*fe2*1;
        fe_mat(stim_curr,:) = p.ae*fe;
        fe2_mat(stim_curr,:) = p.ae2*fe2;
        gi_mat(stim_curr,:)  = p.ai*fi*spfr_data.jumpOn +...
            p.ai2*fi2*spfr_data.jumpOff;
        %gi_mat(stim_curr,:)  = p.ai*fi*1 +...
        %    p.ai2*fi2*1;
        fi_mat(stim_curr,:) = p.ai*fi;
        fi2_mat(stim_curr,:) = p.ai2*fi2;
    end
    
    ge = sum(ge_mat,1);   %sum together all traces over the time course
    gi = sum(gi_mat,1);
    fe = sum(fe_mat,1);
    fe2 = sum(fe2_mat,1);
    fi = sum(fi_mat,1);
    fi2 = sum(fi2_mat,1);
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    %{
      V = V(1:size(spfr_data.time,1));
      ge = ge(1:size(spfr_data.time,1));
      gi = gi(1:size(spfr_data.time,1));
    %}
       
    function f = fc_sol(t_tot, t_on, t_off, Tr, Td)
        t0 = t_tot(t_on);
        t_tot = t_tot - t0; %shift everything so t0 is 0
        t1 = t_tot(t_off);  %so now t1 is effDur
        
        %solve where t0 = 0, so shift everything over
        F1 =...
            @(t) (Td - Tr + Tr.*exp(-t./Tr))./(Td - Tr) -...
            (Td.*exp(-t./Td))./(Td - Tr);
        F2 =...
            @(t) exp(t1./Td).*exp(-t./Td).*...
            ((Td - Tr + Tr.*exp(-t1./Tr))./(Td - Tr) -...
            (Td.*exp(-t1./Td))./(Td - Tr) +...
            (Tr.*exp(t1./Td - t1./Tr).*...
            exp(-t1./Td).*(exp(t1./Tr) - 1))./(Td - Tr)) -...
            (Tr.*exp(t./Td - t./Tr).*...
            exp(-t./Td).*(exp(t1./Tr) - 1))./(Td - Tr);
        f =...
            [zeros(size(t_tot(1:t_on-1)));F1(t_tot(t_on:t_off-1));...
            F2(t_tot(t_off:end))]';
        
         if any(isnan(f) | isinf(f))
            %if nans, solve the long way. t0 is 0, and t1 is effDur
                H1 = @(t) 1 - exp(-t./Tr);
                %F1 defined above
                
                B = H1(t1);
                C = F1(t1);

                t_tot = t_tot-t1; %shift frame where t1 is 0
                t1 = 0;
                H2 = @(t) B.*exp(-t./Tr);
                F2 = @(t) exp(-t./Td).*(C + (B.*Tr)./(Td - Tr)) - (B.*Tr.*exp(-t./Tr))./(Td - Tr);

            counter = 0;
            t_true = 0;
            while any(isnan(f) | isinf(f)) %if still nans, evaluate each ODE
                counter = counter+1;
                if counter> 100 %no unchecked while loops in cluster
                    break
                end            
                tmp_idx = find(isnan(f) | isinf(f),1) -1;  %find where function NaNs

                t_true = t_true + t1; %make t1 now relative to the old t1, so it is the time of integration until NaN
                t1 = t_tot(tmp_idx) - t_true; %t1 is integration time
                B = H2(t1);
                C = F2(t1);

                H2 = @(t) B.*exp(-t./Tr);
                F2 = @(t) exp(-t./Td).*(C + (B.*Tr)./(Td - Tr)) - (B.*Tr.*exp(-t./Tr))./(Td - Tr);
                f(tmp_idx+1:end) = F2(t_tot(tmp_idx+1:end) - t_tot(tmp_idx)); %add to output vector

            end %repeat while we have nans

        end
    end

    function f = fd_sol(t_tot, t_on, t_off, Ti, Tr, Td, A)
        t0 = t_tot(t_on);
        t_tot = t_tot - t0;
        t1 = t_tot(t_off);
        t_tot = t_tot - t1; %shift frame again so that t1 is 0
        %t1 = 0; %S: this should not affect anything
        
        F2 = @(t) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./...
            (Td - Ti) - (A.*Td.*Ti.*Tr.*exp(t./Td - t./Tr))./...
            (Td - Tr)))./(Td.*Ti - Td.*Tr) -...
            (A.*Td.*Ti.*exp(-t./Td))./(Td.*Ti + Td.*Tr - Ti.*Tr - Td^2);
        f = [zeros(size(t_tot(1:t_off-1)));F2(t_tot(t_off:end))]';
        
        if any(isnan(f) | isinf(f))
            %if nans define all functions explicitly and solve piecewise
            A = b;
            B = 0;
            C = 0;
            I2 = @(t) A.*exp(-t./Ti);
            H2 = @(t) exp(-t./Tr).*(B - (A.*Ti)./(Ti - Tr)) + (A.*Ti.*exp(-t./Ti))./(Ti - Tr);
            F2 = @(t) exp(-t./Td).*(C - ((Td.*Tr.*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr) - (A.*Td.*Ti^2)./(Td - Ti))./(Td.*Ti - Td.*Tr)) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr.*exp(t./Td - t./Tr).*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);
            f = [zeros(size(t_tot(1:t_off-1)));F2(t_tot(t_off:end))]';
            
            counter = 0;
            t_true = t1;
            while any(isnan(f) | isinf(f))
                counter = counter+1;
                if counter> 100 %no unchecked while loops in cluster
                    break
                end            
                tmp_idx = find(isnan(f) | isinf(f),1) -1;  %find where function NaNs
                t_true = t_true+t1; %keep running index of last stopping point
                t1 = t_tot(tmp_idx) - t_true; %t1 is integration time, so difference from current stopping point and last stopping point

                A = I2(t1); %evaluate funciton just before nans. these are new init conds to decay
                B = H2(t1);
                C = F2(t1);

                I2 = @(t) A.*exp(-t./Ti); %these functions evaluate from 0
                H2 = @(t) exp(-t./Tr).*(B - (A.*Ti)./(Ti - Tr)) + (A.*Ti.*exp(-t./Ti))./(Ti - Tr);
                F2 = @(t) exp(-t./Td).*(C - ((Td.*Tr.*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr) - (A.*Td.*Ti^2)./(Td - Ti))./(Td.*Ti - Td.*Tr)) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr.*exp(t./Td - t./Tr).*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);

                f(tmp_idx+1:end) = F2(t_tot(tmp_idx+1:end) - t_tot(tmp_idx)); %add to output vector
            end %repeat while we have nans

        end
    end

end
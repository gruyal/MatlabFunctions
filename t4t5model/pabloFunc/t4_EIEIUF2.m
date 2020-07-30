function [V] = t4_EIEIUF(params,spfr_data)

if spfr_data.val == 1                   %all params stored in "rows" of 4, 
    mu_p  = params(1:4);
    sig_p = params(5:8);
    amp_p = params(9:12);
    tr_p  = params(13:16)*10;
    td_p  = params(17:20)*10;
    be_p  = params(23:4:40);
    bi_p  = params(24:4:40);
    ti_p  = params(41:42)*10;
else
    mu_p  = params([3:4,1:2]); %if NC, flip order of all param sets (seconds always fed into delta)
    sig_p = params([7:8,5:6]);
    amp_p = params([11:12,9:10]);
    tr_p  = params([15:16,13:14])*10;
    td_p  = params([19:20,17:18])*10;
    be_p  = params(21:4:40);
    bi_p  = params(22:4:40);
    ti_p  = params([42,41])*10;
end


p.Vr = 0;
p.Ve = 65;
p.Vi = -10;

%stim parameters
width = spfr_data.width-1;
p.width = spfr_data.width;
p.stimDur = spfr_data.stimDur;
pos_vect = spfr_data.pos_vect;
d = 30/spfr_data.fr; %set manual delay to replicate delay in transmission from photoreceptors to T4/5 response
tmp = find(spfr_data.stimIdx);
stimIdx = false(size(spfr_data.stimIdx));
stimIdx(tmp+d) = true;
stim_time = fix(spfr_data.time(stimIdx));

buffer = 5;
M = 1:(length(spfr_data.pos_vect)+2*buffer);
x = M+min(spfr_data.pos_vect)-buffer;

Ie = zeros(length(stim_time),length(M));
Ii = zeros(length(stim_time),length(M));

for i = 1:length(pos_vect)
    pos = pos_vect(i) - min(pos_vect)+buffer;

     Ie(i,pos-width:pos) = 1;
     Ii(i,pos-width:pos) = 1;

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
if ~spfr_data.cat_flag
    effDur = effDur*spfr_data.width; %if moving, then also account for width 
end

% p.be2 = (1 - exp(-effDur/tt_p(1))); %expo  plateau, N(t) = No*(1- exp(-t/tau)).
% p.bi2 = (1 - exp(-effDur/tt_p(2))); %expo  plateau, N(t) = No*(1- exp(-t/tau)).
switch effDur
    case 40
        p.be2 = be_p(1);
        p.bi2 = bi_p(1);
    case 80
        p.be2 = be_p(2);
        p.bi2 = bi_p(2);
    case 160
        p.be2 = be_p(3);
        p.bi2 = bi_p(3);
    case 320
        p.be2 = be_p(4);
        p.bi2 = bi_p(4);
    case 640
        p.be2 = be_p(5);
        p.bi2 = bi_p(5);        
end
p.Tic = ti_p(1);
p.Tid = ti_p(2);


% Model
p.ae = Ae.*exp(-(x-mue).^2 / (2*sige^2) );
p.ai = Ai.*exp(-(x-mui).^2 / (2*sigi^2) );
p.ae2 = Ae2.*exp(-(x-mue2).^2 / (2*sige2^2) );
p.ai2 = Ai2.*exp(-(x-mui2).^2 / (2*sigi2^2) );

p.stim_time = stim_time;
y0 = zeros((8*size(Ie,2)),1);
t_ind = find([stimIdx;1]);
[~,s_ind] = min(abs(spfr_data.time - (spfr_data.time(stimIdx) + effDur)'));
s_ind = s_ind';

p.stim_time = stim_time;

    %populate a matrix with this dynamic for each stim
    fe = zeros(size(p.Ie,2),size(spfr_data.time,1));
    fi = zeros(size(p.Ii,2),size(spfr_data.time,1));
    fe2 = zeros(size(p.Ii,2),size(spfr_data.time,1));
    fi2 = zeros(size(p.Ii,2),size(spfr_data.time,1));

if spfr_data.cat_flag     
    %subtract delay from ending index to not encroach on subsequent SPFR
    t_tot = spfr_data.time(1:t_ind(2)-1);
    t_on = t_ind(1);
    t_off= s_ind(1);
    
    fe_tmp = fc_sol(t_tot, t_on, t_off, p.Tic, p.Tre, p.Tde);
    fi_tmp = fc_sol(t_tot, t_on, t_off, p.Tic, p.Tri, p.Tdi);
    
    fe2_tmp = fd_sol(t_tot, t_on, t_off, p.Tid, p.Tre2, p.Tde2, p.be2);
    fi2_tmp = fd_sol(t_tot, t_on, t_off, p.Tid, p.Tri2, p.Tdi2, p.bi2);
   
    
    %after solving for conductances at each time, we can now shift back
    %to non-delayed frame to fill vector
    t_ind = t_ind-d;
    s_ind = s_ind-d;
    for i = 1:size(Ie,1)
        fe(logical(p.Ie(i,:)),t_ind(i):(t_ind(i)+length(fe_tmp)-1)) = repmat(fe_tmp,sum(p.Ie(i,:)),1);
        fi(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+length(fi_tmp)-1)) = repmat(fi_tmp,sum(p.Ii(i,:)),1);
        fe2(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+length(fe2_tmp)-1)) = repmat(fe2_tmp,sum(p.Ii(i,:)),1);
        fi2(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+length(fe2_tmp)-1)) = repmat(fi2_tmp,sum(p.Ii(i,:)),1);
    end
    
    %caluclate conductances and steady state voltage
    ge  = p.ae*fe + p.ae2*fe2;
    gi  = p.ai*fi + p.ai2*fi2;
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = V(1:size(spfr_data.time,1));
    
else
    %subtract delay from ending index to not encroach on subsequent SPFR
    %spfr_data.time = spfr_data.time - spfr_data.time(t_ind(1)); %because our analytical solutions assume t0 = 0
    t_tot = spfr_data.time;
    
    t_ind = [t_ind(1)*ones(width,1);t_ind];
    effDur = spfr_data.stimDur;
    effDur = effDur*sum(p.Ie,1); %if moving, then also account for width 
    fe_tmp = zeros(length(unique(effDur)),length(t_tot));
    fi_tmp = zeros(length(unique(effDur)),length(t_tot));
    fe2_tmp = zeros(length(unique(effDur)),length(t_tot));
    fi2_tmp = zeros(length(unique(effDur)),length(t_tot));

    for effDuri = unique(effDur)                        
                switch effDuri
                    case 0 
                        continue
                    case 40
                        be2 = be_p(1);
                        bi2 = bi_p(1);
                    case 80
                        be2 = be_p(2);
                        bi2 = bi_p(2);
                    case 160
                        be2 = be_p(3);
                        bi2 = bi_p(3);
                    case 320
                        be2 = be_p(4);
                        bi2 = bi_p(4);
                    case 640
                        be2 = be_p(5);
                        bi2 = bi_p(5);        
                end
            
            [~,s] = min(abs(spfr_data.time - spfr_data.time(t_ind(1)) - effDuri));
            t_on = t_ind(1);
            t_off = s;
            fe_tmp(1 + effDuri/p.stimDur,:) = fc_sol(t_tot, t_on, t_off, p.Tic, p.Tre, p.Tde);
            fi_tmp(1 + effDuri/p.stimDur,:) = fc_sol(t_tot, t_on, t_off, p.Tic, p.Tri, p.Tdi);
            fe2_tmp(1 + effDuri/p.stimDur,:) = fd_sol(t_tot, t_on, t_off, p.Tid, p.Tre2, p.Tde2, be2);
            fi2_tmp(1 + effDuri/p.stimDur,:) = fd_sol(t_tot, t_on, t_off, p.Tid, p.Tri2, p.Tdi2, bi2);        
    end
    
    
    t_ind = t_ind - t_ind(1)+1;
    ii = 0;
    for i = 1:size(p.Ie,2)
        if ~any(p.Ie(:,i),1)
            continue
        end
        ii = ii+1;  
        fe(i,t_ind(ii):end)  =  fe_tmp(1+effDur(i)/p.stimDur, 1:end-t_ind(ii)+1);
        fi(i,t_ind(ii):end)  =  fi_tmp(1+effDur(i)/p.stimDur, 1:end-t_ind(ii)+1);
        fe2(i,t_ind(ii):end) = fe2_tmp(1+effDur(i)/p.stimDur, 1:end-t_ind(ii)+1);
        fi2(i,t_ind(ii):end) = fi2_tmp(1+effDur(i)/p.stimDur, 1:end-t_ind(ii)+1);
    end
    
    if mean(diff(pos_vect)) < 0 %if ND, we need to reverse order
        fe(any(fe,2),:) = flipud(fe(any(fe,2),:));
        fi(any(fi,2),:) = flipud(fi(any(fi,2),:));
        fe2(any(fe2,2),:) = flipud(fe2(any(fe2,2),:));
        fi2(any(fi2,2),:) = flipud(fi2(any(fi2,2),:));
    end
    
    ge  = p.ae*fe + p.ae2*fe2;
    gi  = p.ai*fi + p.ai2*fi2;  
    
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = [zeros(length(spfr_data.baseSub) - length(V),1);V]; %because of delay, add bad time with zeros

end

    function f = fc_sol(t_tot, t_on, t_off, Ti, Tr, Td)

        t0 = t_tot(t_on);
        t_tot = t_tot - t0; %shift everything so t0 is 0
        t1 = t_tot(t_off);  %so now t1 is effDur

        %solve where t0 = 0, so shift everything over
        F1 = @(t) (exp(-t./Td).*(Td.*exp(t./Td).*(Ti - Tr) + (Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr^2.*exp(t./Td - t./Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr) + (Td^2.*exp(-t./Td))./(Td.*Ti + Td.*Tr - Ti.*Tr - Td^2);
        F2 = @(t) (exp(-t./Td).*(Td^2 - Td^2.*exp(t1./Td)))./(Td.*Ti + Td.*Tr - Ti.*Tr - Td^2) + (exp(-t./Td).*((Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr^2.*exp(t./Td - t./Tr))./(Td - Tr) - (Td.*Ti^2.*exp(t./Td - t./Ti + t1./Ti))./(Td - Ti) + (Td.*Tr^2.*exp(t./Td - t./Tr + t1./Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);
        f = [zeros(size(t_tot(1:t_on-1)));F1(t_tot(t_on:t_off-1));F2(t_tot(t_off:end))]';

        if any(isnan(f) | isinf(f))
            %if nans, solve the long way. t0 is 0, and t1 is effDur
                I1 = @(t) 1 - exp(-t/Ti);
                H1 = @(t) (Tr*exp(-t/Tr))/(Ti - Tr) - (Tr - Ti + Ti*exp(-t/Ti))/(Ti - Tr);
                F1 = @(t) (exp(-t/Td)*(Td*exp(t/Td)*(Ti - Tr) + (Td*Ti^2*exp(t/Td - t/Ti))/(Td - Ti) - (Td*Tr^2*exp(t/Td - t/Tr))/(Td - Tr)))/(Td*Ti - Td*Tr) + (Td^2*exp(-t/Td))/(Td*Ti + Td*Tr - Ti*Tr - Td^2);

                A = I1(t1); %find init conds for decay
                B = H1(t1);
                C = F1(t1);

                t_tot = t_tot-t1; %shift frame where t1 is 0
                t1 = 0;
                I2 = @(t) A.*exp(-t./Ti); %evaluate decays from this time (t1 = 0)
                H2 = @(t) exp(-t./Tr).*(B - (A.*Ti)./(Ti - Tr)) + (A.*Ti.*exp(-t./Ti))./(Ti - Tr);
                F2 = @(t) exp(-t./Td).*(C - ((Td.*Tr.*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr) - (A.*Td.*Ti^2)./(Td - Ti))./(Td.*Ti - Td.*Tr)) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr.*exp(t./Td - t./Tr).*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);

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
                A = I2(t1); %find initial conditions at this time
                B = H2(t1);
                C = F2(t1);

                I2 = @(t) A.*exp(-t./Ti); %evaluate decays from this time
                H2 = @(t) exp(-t./Tr).*(B - (A.*Ti)./(Ti - Tr)) + (A.*Ti.*exp(-t./Ti))./(Ti - Tr);
                F2 = @(t) exp(-t./Td).*(C - ((Td.*Tr.*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr) - (A.*Td.*Ti^2)./(Td - Ti))./(Td.*Ti - Td.*Tr)) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr.*exp(t./Td - t./Tr).*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);
                f(tmp_idx+1:end) = F2(t_tot(tmp_idx+1:end) - t_tot(tmp_idx)); %add to output vector

            end %repeat while we have nans

        end
    end

    function f = fd_sol(t_tot, t_on, t_off, Ti, Tr, Td, b)

        t0 = t_tot(t_on);
        t_tot = t_tot - t0;
        t1 = t_tot(t_off);

        A = b;

        t_tot = t_tot - t1; %shift frame again so that t1 is 0
        t1 = 0;

        F2 = @(t) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (A.*Td.*Ti.*Tr.*exp(t./Td - t./Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr) - (A.*Td.*Ti.*exp(-t./Td))./(Td.*Ti + Td.*Tr - Ti.*Tr - Td^2);
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
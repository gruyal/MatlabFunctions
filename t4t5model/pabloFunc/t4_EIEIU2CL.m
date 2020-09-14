function [V,ge,gi,g_mat,s_mat,num,den] = t4_EIEIU2CL(params,spfr_data)

if spfr_data.val == 1                   %all params stored in "rows" of 4, 
    mu_p  = params(1:4);
    sig_p = params(5:8);
    amp_p = params(9:12);
    tr_p  = params([13:16])*10;
    td_p  = params([17:20])*10;
    be_p  = params(23)*.01;
    bi_p  = params(24)*.01;
    ti_p  = params(26)*10;
else
    mu_p  = params([3:4,1:2]); %if NC, flip order of all param sets (seconds always fed into delta)
    sig_p = params([7:8,5:6]);
    amp_p = params([11:12,9:10]);
    tr_p  = params([15:16,13:14])*10;
    td_p  = params([19:20,17:18])*10;
    be_p  = params(21)*.01;
    bi_p  = params(22)*.01;
    ti_p  = params(25)*10;
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

%regardless of contrast, both linear
p.be2 = effDur*be_p;
p.bi2 = effDur*bi_p;


p.Tid = ti_p;

% Model
p.ae = Ae.*exp(-(x-mue).^2 / (2*sige^2) );
p.ai = Ai.*exp(-(x-mui).^2 / (2*sigi^2) );
p.ae2 = Ae2.*exp(-(x-mue2).^2 / (2*sige2^2) );
p.ai2 = Ai2.*exp(-(x-mui2).^2 / (2*sigi2^2) );

p.stim_time = stim_time;
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
    
    fe_tmp = fc_sol(t_tot, t_on, t_off, p.Tre, p.Tde);
    fi_tmp = fc_sol(t_tot, t_on, t_off, p.Tri, p.Tdi);
    
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
    
    fe = fe(:,1:size(spfr_data.time,1));
    fi = fi(:,1:size(spfr_data.time,1));
    fe2 = fe2(:,1:size(spfr_data.time,1));
    fi2 = fi2(:,1:size(spfr_data.time,1));
    
    %caluclate conductances and steady state voltage
    ge  = p.ae*fe + p.ae2*fe2;
    gi  = p.ai*fi + p.ai2*fi2;
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = V(1:size(spfr_data.time,1));
    num = (p.Vr + ge*p.Ve + gi*p.Vi);
    den = (1+ge+gi);
    ge1 = (p.ae*fe).*(65-V)';
    gi1 = (p.ai*fi).*(10+V)';
    ge2 = (p.ae2*fe2).*(65-V)';
    gi2 = (p.ai2*fi2).*(10+V)';
    if spfr_data.val == 0
        g_mat = {ge1,gi1,ge2,gi2};
        s_mat = {mue,mui,mue2,mui2;...
                 sige,sigi,sige2,sigi2};
    else
        g_mat = {ge2,gi2,ge1,gi1};
        s_mat = {mue2,mui2,mue,mui;...
                 sige2,sigi2,sige,sigi};
    end
    
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

    counter2 = 0;
    for effDuri = unique(effDur)
            counter2 = counter2+1;
                        
            be2 = effDuri*be_p; %if PC, deltas are constants.
            bi2 = effDuri*bi_p; 
                  
            [~,s] = min(abs(spfr_data.time - spfr_data.time(t_ind(1)) - effDuri));
            t_on = t_ind(1);
            t_off = s;
            fe_tmp(1 + effDuri/p.stimDur,:) = fc_sol(t_tot, t_on, t_off, p.Tre, p.Tde);
            fi_tmp(1 + effDuri/p.stimDur,:) = fc_sol(t_tot, t_on, t_off, p.Tri, p.Tdi);
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
    
    fe = fe(:,1:size(spfr_data.time,1));
    fi = fi(:,1:size(spfr_data.time,1));
    fe2 = fe2(:,1:size(spfr_data.time,1));
    fi2 = fi2(:,1:size(spfr_data.time,1));
    
    %caluclate conductances and steady state voltage
   
    ge  = p.ae*fe + p.ae2*fe2;
    gi  = p.ai*fi + p.ai2*fi2;  
    
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = [zeros(length(spfr_data.baseSub) - length(V),1);V]; %because of delay, add bad time with zeros

    num = (p.Vr + ge*p.Ve + gi*p.Vi);
    den = (1+ge+gi);
    ge1 = (p.ae*fe).*(65-V)';
    gi1 = (p.ai*fi).*(10+V)';
    ge2 = (p.ae2*fe2).*(65-V)';
    gi2 = (p.ai2*fi2).*(10+V)';
    if spfr_data.val == 1
        g_mat = {ge1,gi1,ge2,gi2};
        s_mat = {mue,mui,mue2,mui2;...
                 sige,sigi,sige2,sigi2};
    else
        g_mat = {ge2,gi2,ge1,gi1};
        s_mat = {mue2,mui2,mue,mui;...
                 sige2,sigi2,sige,sigi};
    end
end


    function f = fc_sol(t_tot, t_on, t_off, Tr, Td)
        if Tr == Td
            Tr = Tr + 1e-7; %add a step size smaller than step tolerance of optimizer but avoids discontinuity
        end
        
        t0 = t_tot(t_on);
        t_tot = t_tot - t0; %shift everything so t0 is 0
        t1 = t_tot(t_off);  %so now t1 is effDur

        %solve where t0 = 0, so shift everything over
        F1 = @(t) (Td - Tr + Tr.*exp(-t./Tr))./(Td - Tr) - (Td.*exp(-t./Td))./(Td - Tr);
        F2 = @(t) exp(t1./Td).*exp(-t./Td).*((Td - Tr + Tr.*exp(-t1./Tr))./(Td - Tr) - (Td.*exp(-t1./Td))./(Td - Tr) + (Tr.*exp(t1./Td - t1./Tr).*exp(-t1./Td).*(exp(t1./Tr) - 1))./(Td - Tr)) - (Tr.*exp(t./Td - t./Tr).*exp(-t./Td).*(exp(t1./Tr) - 1))./(Td - Tr);
%         if Tr == Td
%             F1 = @(t) (Tr.*exp(-t./Tr).*(2.*Ti - Tr))./(Ti - Tr)^2 + (exp(-t./Tr).*(Tr.*t + Tr.*exp(t./Tr).*(Ti - Tr) - (Ti^2.*Tr.*exp((t.*(Ti - Tr))./(Ti.*Tr)))./(Ti - Tr)))./(Tr.*(Ti - Tr));
%             F2 = @(t) (exp(-t./Tr).*(Tr^2.*exp(t1./Tr) + 2.*Ti.*Tr - Tr^2 + Ti.*t1.*exp(t1./Tr) - Tr.*t1.*exp(t1./Tr) - 2.*Ti.*Tr.*exp(t1./Tr)))./(Ti - Tr)^2 + (exp(-t./Tr).*(Tr.*t - Tr.*t.*exp(t1./Tr) + (Ti^2.*Tr.*exp((Ti.*t - Tr.*t + Tr.*t1)./(Ti.*Tr)))./(Ti - Tr) - (Ti^2.*Tr.*exp((t.*(Ti - Tr))./(Ti.*Tr)))./(Ti - Tr)))./(Tr.*(Ti - Tr));
%         end       
        
        f = [zeros(size(t_tot(1:t_on-1)));F1(t_tot(t_on:t_off-1));F2(t_tot(t_off:end))]';

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

    function f = fd_sol(t_tot, t_on, t_off, Ti, Tr, Td, b)
        if Tr == Td
            Tr = Tr + 1e-7; %add a step size smaller than step tolerance of optimizer but avoids discontinuity
        end
        t0 = t_tot(t_on);
        t_tot = t_tot - t0;
        t1 = t_tot(t_off);

        A = b;

        t_tot = t_tot - t1; %shift frame again so that t1 is 0
        t1 = 0;

        F2 = @(t) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (A.*Td.*Ti.*Tr.*exp(t./Td - t./Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr) - (A.*Td.*Ti.*exp(-t./Td))./(Td.*Ti + Td.*Tr - Ti.*Tr - Td^2);
%         if Tr == Td
%            F2 = @(t) -(A.*Ti.*exp(-t./Tr).*(Ti.*Tr + Ti.*t - Tr.*t - Ti.*Tr.*exp((t.*(Ti - Tr))./(Ti.*Tr))))./(Tr.*(Ti - Tr)^2);
%         end
        
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
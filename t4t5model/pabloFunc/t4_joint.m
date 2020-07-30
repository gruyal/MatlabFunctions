function [V] = t4_EIEI_joint(params,spfr_data)

if spfr_data.val == 0                   %all params stored in "rows" of 4, 
    mu_p  = params(1:4);
    sig_p = params(5:8);
    amp_p = params(9:12);
    tr_p  = params(17:20)*10;
    td_p  = params(21:24)*10;
    tt_p  = params(33:34)*10;
else
    mu_p  = params([3:4,1:2]); %if NC, flip order of all param sets (seconds always fed into delta)
    sig_p = params([7:8,5:6]);
    amp_p = params([15:16,13:14]);
    tr_p  = params([27:28,25:26])*10;
    td_p  = params([31:32,29:30])*10;
    tt_p  = params(35:36)*10; %and take tau for second delta
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
spfr_data.stimIdx = false(size(spfr_data.stimIdx));
spfr_data.stimIdx(tmp+d) = true;
stim_time = fix(spfr_data.time(spfr_data.stimIdx));

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

p.be2 = (1 - exp(-effDur/tt_p(1))); %expo  plateau, N(t) = No*(1- exp(-t/tau)).
p.bi2 = (1 - exp(-effDur/tt_p(2))); %expo  plateau, N(t) = No*(1- exp(-t/tau)).



% Model
p.ae = Ae.*exp(-(x-mue).^2 / (2*sige^2) );
p.ai = Ai.*exp(-(x-mui).^2 / (2*sigi^2) );
p.ae2 = Ae2.*exp(-(x-mue2).^2 / (2*sige2^2) );
p.ai2 = Ai2.*exp(-(x-mui2).^2 / (2*sigi2^2) );

p.stim_time = stim_time;
y0 = zeros((8*size(Ie,2)),1);
t_ind = find([spfr_data.stimIdx;1]);
[~,s_ind] = min(abs(spfr_data.time - (spfr_data.time(spfr_data.stimIdx) + p.stimDur)'));
s_ind = s_ind';

p.stim_time = stim_time;

if spfr_data.cat_flag   
    
    %excitatory (classic)
    he_fun1 = @(t) 1-exp(-t/p.Tre);

    if abs(p.Tre - p.Tde) >= 1e8*(p.Tre - p.Tde)^2
        global taylor_flag
        taylor_flag = true;
        E = p.Tre - p.Tde;
        fe_fun1 = @(t) 1 - exp(-t/p.Tde) - (1/p.Tde).*exp(-t/p.Tde).*(t + (E*t.^2)./(2*p.Tde^2));
        fe_fun2 = @(t) exp(-t./p.Tde).*(fe_fun1(p.stimDur) + (he_fun1(p.stimDur)./p.Tde).*(t - E*(t+p.stimDur).*(t)./(p.Tde^2)));
    else
        a = (p.Tre - p.Tde)/(p.Tre*p.Tde);
        fe_fun1 = @(t) (exp(-t/p.Tde)/p.Tde).*(p.Tde*exp(t/p.Tde) - p.Tde - (1/a)*exp(a*t) + (1/a));
        fe_fun2 = @(t) exp(-t/p.Tde).*(fe_fun1(p.stimDur) + he_fun1(p.stimDur)/(a*p.Tde)*(exp(a*t) -1));
        if any(isnan(fe_fun2(spfr_data.time(s_ind(1):t_ind(2)-1)))) || any(isinf(fe_fun2(spfr_data.time(s_ind(1):t_ind(2)-1))))
            fe_fun1 = @(t) (exp(-t/p.Tde)/p.Tde).*(p.Tde*exp(t/p.Tde) - p.Tde - (1/a)*exp(a*t) + (1/a));
            fe_fun2 = @(t) exp(-t/p.Tde+500).*(fe_fun1(p.stimDur)*exp(-500) + he_fun1(p.stimDur)/(a*p.Tde)*(exp((a*t - 500)) -exp(-500)));
        end
    end   
    
    %inhibitory (classic)
    hi_fun1 = @(t) 1-exp(-t/p.Tri);
 
    if abs(p.Tri - p.Tdi) >= 1e8*(p.Tri - p.Tdi)^2
        global taylor_flag
        taylor_flag = true;
        E = p.Tri - p.Tdi;
        fi_fun1 = @(t) 1 - exp(-t/p.Tdi) - (1/p.Tdi).*exp(-t/p.Tdi).*(t + (E*t.^2)./(2*p.Tdi^2));
        fi_fun2 = @(t) exp(-t./p.Tdi).*(fi_fun1(p.stimDur) + (hi_fun1(p.stimDur)./p.Tdi).*(t - E*(t+p.stimDur).*(t)./(p.Tdi^2)));
    else
        a = (p.Tri - p.Tdi)/(p.Tri*p.Tdi);
        fi_fun1 = @(t) (exp(-t/p.Tdi)/p.Tdi).*(p.Tdi*exp(t/p.Tdi) - p.Tdi - (1/a)*exp(a*t) + (1/a));
        fi_fun2 = @(t) exp(-t/p.Tdi).*(fi_fun1(p.stimDur) + hi_fun1(p.stimDur)/(a*p.Tdi)*(exp(a*t) -1));
        if any(isnan(fi_fun2(spfr_data.time(s_ind(1):t_ind(2)-1)))) || any(isinf(fi_fun2(spfr_data.time(s_ind(1):t_ind(2)-1))))
            fi_fun1 = @(t) (exp(-t/p.Tdi)/p.Tdi).*(p.Tdi*exp(t/p.Tdi) - p.Tdi - (1/a)*exp(a*t) + (1/a));
            fi_fun2 = @(t) exp(-t/p.Tdi+500).*(fi_fun1(p.stimDur)*exp(-500) + hi_fun1(p.stimDur)/(a*p.Tdi)*(exp((a*t - 500)) -exp(-500)));
        end
    end
    
    %excitatory2 (delta)
    if abs(p.Tre2 - p.Tde2) >= 1e8*(p.Tre2 - p.Tde2)^2
        global taylor_flag
        taylor_flag = true;
        E = p.Tre2 - p.Tde2;
        fe2_fun2 = @(t) exp(-t./p.Tde2).*(p.be2./p.Tde2).*(t - E*(t+p.stimDur).*(t)./(p.Tde2^2));
    else
        a = (p.Tre2 - p.Tde2)/(p.Tre2*p.Tde2);
        fe2_fun2 = @(t) exp(-t/p.Tde2).*(p.be2/(p.Tde2*a)).*(exp(a*t)-1);
        if any(isnan(fe2_fun2(spfr_data.time(s_ind(1):t_ind(2)-1)))) || any(isinf(fe2_fun2(spfr_data.time(s_ind(1):t_ind(2)-1))))
            fe2_fun2 = @(t) exp(-t/p.Tde2+500).*((p.be2/(a*p.Tde2))*(exp((a*t - 500)) -exp(-500)));
        end  
    end
    
    %inhibitory2 (delta)
    if abs(p.Tri2 - p.Tdi2) >= 1e8*(p.Tri2 - p.Tdi2)^2
        global taylor_flag
        taylor_flag = true;
        E = p.Tri2 - p.Tdi2;
        fi2_fun2 = @(t) exp(-t./p.Tdi2).*(p.bi2./p.Tdi2).*(t - E*(t+p.stimDur).*(t)./(p.Tdi2^2));
    else
        a = (p.Tri2 - p.Tdi2)/(p.Tri2*p.Tdi2);
        fi2_fun2 = @(t) exp(-t/p.Tdi2).*(p.bi2/(p.Tdi2*a)).*(exp(a*t)-1);
        if any(isnan(fi2_fun2(spfr_data.time(s_ind(1):t_ind(2)-1)))) || any(isinf(fi2_fun2(spfr_data.time(s_ind(1):t_ind(2)-1))))
            fi2_fun2 = @(t) exp(-t/p.Tdi2+500).*((p.bi2/(a*p.Tdi2))*(exp((a*t - 500)) -exp(-500)));
        end  
    end
    
    %subtract delay from ending index to not encroach on subsequent SPFR
    spfr_data.time = spfr_data.time - spfr_data.time(t_ind(1)); %because our analytical solutions assume t0 = 0

    fe_tmp = zeros(size(spfr_data.time(t_ind(1):t_ind(2)-1)));
    fe_tmp(t_ind(1):s_ind(1)-1) = fe_fun1(spfr_data.time(t_ind(1):s_ind(1)-1));
    fe_tmp(s_ind(1):t_ind(2)-1-d) = fe_fun2(spfr_data.time(s_ind(1):t_ind(2)-1-d) - p.stimDur);
    
    fi_tmp = zeros(size(spfr_data.time(t_ind(1):t_ind(2)-1)));
    fi_tmp(t_ind(1):s_ind(1)-1) = fi_fun1(spfr_data.time(t_ind(1):s_ind(1)-1));
    fi_tmp(s_ind(1):t_ind(2)-1-d) = fi_fun2(spfr_data.time(s_ind(1):t_ind(2)-1-d) - p.stimDur);
    
    fe2_tmp = zeros(size(spfr_data.time(t_ind(1):t_ind(2)-1)));
    fe2_tmp(s_ind(1):t_ind(2)-1-d) = fe2_fun2(spfr_data.time(s_ind(1):t_ind(2)-1-d) - p.stimDur);
    
    fi2_tmp = zeros(size(spfr_data.time(t_ind(1):t_ind(2)-1)));
    fi2_tmp(s_ind(1):t_ind(2)-1-d) = fi2_fun2(spfr_data.time(s_ind(1):t_ind(2)-1-d) - p.stimDur);
    
    %populate a matrix with this dynamic for each stim
    fe = zeros(size(p.Ie,2),size(spfr_data.time,1));
    fi = zeros(size(p.Ii,2),size(spfr_data.time,1));
    fe2 = zeros(size(p.Ii,2),size(spfr_data.time,1));
    fi2 = zeros(size(p.Ii,2),size(spfr_data.time,1));
    
    %after solving for conductances at each time, we can now shift back
    %to non-delayed frame to fill vector
    t_ind = t_ind-d;
    s_ind = s_ind-d;
    for i = 1:size(Ie,1)
        fe(logical(p.Ie(i,:)),t_ind(i):(t_ind(i)+size(fe_tmp,1)-1)) = repmat(fe_tmp,1,sum(p.Ie(i,:)))';
        fi(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+size(fi_tmp,1)-1)) = repmat(fi_tmp,1,sum(p.Ii(i,:)))';
        fe2(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+size(fe2_tmp,1)-1)) = repmat(fe2_tmp,1,sum(p.Ii(i,:)))';
        fi2(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+size(fe2_tmp,1)-1)) = repmat(fi2_tmp,1,sum(p.Ii(i,:)))';
    end
    
    %caluclate conductances and steady state voltage
    ge  = p.ae*fe + p.ae2*fe2;
    gi  = p.ai*fi + p.ai2*fi2;
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = V(1:size(spfr_data.time,1));
    
    
    if any(isnan(V))
        global nan_flag
        nan_flag = true;
    end
    
else
    if any(params(17:end)*10 < .5)
        V = zeros(size(spfr_data.baseSub));
    else
    
    p.Ie = [p.Ie;zeros(2,size(p.Ie,2))];
    t_ind = [t_ind(1:end-1);s_ind(end);t_ind(end)];
    tmp = cell(length(t_ind)-1,1);
    effDur = spfr_data.stimDur;

    for i = 1:length(t_ind)-1
        sol = ode45(@(t,y) vm_wrap_mov(t,y,p),[spfr_data.time(t_ind(i)) spfr_data.time(t_ind(i+1)-1)+10],y0);
        tspan = spfr_data.time(t_ind(i):t_ind(i+1)-1); % solve during bar flash
        tmp{i} = deval(sol,tspan); %solve each bit and store laterally
        gate = (p.Ie(i,:) - p.Ie(i+1,:)) > 0;
        y0 = tmp{i}(:,end); %next bit picks up where this left off
        effDur2 = effDur*sum(p.Ie(:,gate));
        p.be2 = (1 - exp(-effDur2/tt_p(1))); %expo  plateau, N(t) = No*(1- exp(-t/tau)).
        p.bi2 = (1 - exp(-effDur2/tt_p(2))); %expo  plateau, N(t) = No*(1- exp(-t/tau)).
        tmpe = zeros(size(gate))';
        tmpi = zeros(size(gate))';
        tmpe(gate) = p.be2;
        tmpi(gate) = p.bi2;
        y0(5:8:length(y0),end) = tmpe + tmp{i}(5:8:length(y0),end).*~gate'; % set delta function, which jumps to amplitude instantly only at section that relaxes at t_ind
        y0(7:8:length(y0),end) = tmpi + tmp{i}(7:8:length(y0),end).*~gate';
    end
    tmp = [tmp{:}];
    fe = tmp(2:8:length(y0),:);
    fi = tmp(4:8:length(y0),:);
    fe2 = tmp(6:8:length(y0),:);
    fi2 = tmp(8:8:length(y0),:);
    
    ge  = p.ae*fe + p.ae2*fe2;
    gi  = p.ai*fi + p.ai2*fi2;  
    
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = [zeros(length(spfr_data.baseSub) - length(V),1);V]; %because of delay, add bad time with zeros

    end
end
end

%ODE45
function dydt = vm_wrap_mov(t,y,p)
    p.stim = 0;
    p.stim_num = 1;
    id1 = t > p.stim_time ;
    id2 = t < p.stim_time + p.stimDur;

    tmp = find(id1.*id2);
    if any(tmp)
        p.stim = 1;
        p.stim_num = tmp(1);
    end
    dydt = vm_solve(y,p);
end


function dydt = vm_solve(y,p)
he = y(1:8:end);
fe = y(2:8:end);
hi = y(3:8:end);
fi = y(4:8:end);
he2 = y(5:8:end);
fe2 = y(6:8:end);
hi2 = y(7:8:end);
fi2 = y(8:8:end);

dhe = (-he' + p.stim*p.Ii(p.stim_num,:)) / p.Tre;
dfe = (-fe' + he')     / p.Tde;
dhi = (-hi' + p.stim*p.Ii(p.stim_num,:)) / p.Tri;
dfi = (-fi' + hi')     / p.Tdi;
dhe2 = (-he2') / p.Tre2;
dfe2 = (-fe2' + he2')     / p.Tde2;
dhi2 = (-hi2') / p.Tri2;
dfi2 = (-fi2' + hi2')     / p.Tdi2;

tmp = [dhe;dfe;dhi;dfi;dhe2;dfe2;dhi2;dfi2];
dydt = reshape(tmp,[],1);

end
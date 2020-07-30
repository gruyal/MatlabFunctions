function [V] = t5_deltaEA(params,spfr_data)

params = params.*[10,10,1,1,1,10,10,1,1,1];

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

M = 1:(length(spfr_data.pos_vect));
x = min(pos_vect):max(pos_vect);

%if moving bar, we need to tag stims to end of pos_vect
if ~spfr_data.cat_flag
    d = mean(diff(pos_vect)) > 0; %if diff >1, then it's PD, if it's less than 1 it's ND
    if d && width > 0
        pos_vect = [pos_vect, pos_vect(end)+(1:(width))]; %append stim positions to PD side of pos_vect
    elseif ~d && width > 0
        pos_vect = [pos_vect(1)+((width):-1:1), pos_vect];

    end
else %if spfr, we dont want a window
    x =  [x(1)+(-width:-1), x];
end
    
Ie = zeros(length(stim_time),length(M));
Ii = zeros(length(stim_time),length(M));

%pos_vect defines leadingedge of stimulated position. x defines valid
%positions in window
for i = 1:length(pos_vect)
    idx = x <= pos_vect(i) & x >= pos_vect(i)-width;

     Ie(i,idx) = 1;
     Ii(i,idx) = 1;

end

p.Ie = Ie;
p.Ii = Ii;

%excitatory
p.Tre = params(1);
p.Tde = params(2);
Ae = 1;
mue = params(4);
sige = params(5);

%inhibitory
p.Tri = params(6); 
p.Tdi = params(7);   
Ai = params(8); 
mui = params(9);   
sigi = params(10);  


%time
effDur = spfr_data.stimDur;
if ~spfr_data.cat_flag
    effDur = effDur*spfr_data.width; %if moving, then also account for width 
end

p.b = params(3)*(1 - exp(-effDur/params(1))); %expo  plateau, N(t) = No*(1- exp(-t/tau)).

% Model
p.ae = Ae.*exp(-(x-mue).^2 / (2*sige^2) );
p.ai = Ai.*exp(-(x-mui).^2 / (2*sigi^2) );

p.stim_time = stim_time;
y0 = zeros((4*size(Ie,2)),1);
t_ind = find([spfr_data.stimIdx;1]);
[~,s_ind] = min(abs(spfr_data.time - (spfr_data.time(spfr_data.stimIdx) + p.stimDur)'));
s_ind = s_ind';

p.stim_time = stim_time;

if spfr_data.cat_flag     
   
    if abs(p.Tre - p.Tde) >= 1e8*(p.Tre - p.Tde)^2
        global taylor_flag
        taylor_flag = true;
        E = p.Tre - p.Tde;
        fe_fun2 = @(t) exp(-t./p.Tdi).*(p.b./p.Tdi).*(t - E*(t+p.stimDur).*(t)./(p.Tdi^2));
    else
        a = (p.Tre - p.Tde)/(p.Tre*p.Tde);
        fe_fun2 = @(t) exp(-t/p.Tde).*(p.b/(p.Tde*a)).*(exp(a*t)-1);
        if any(isnan(fe_fun2(spfr_data.time(s_ind(1):t_ind(2)-1)))) || any(isinf(fe_fun2(spfr_data.time(s_ind(1):t_ind(2)-1))))
            fe_fun2 = @(t) exp(-t/p.Tde+500).*((p.b/(a*p.Tde))*(exp(a*t - 500) -exp(-500)));
        end  
    end
    
    hi_fun1 = @(t) 1-exp(-t/p.Tri);
    %hi_fun2 = @(t) exp(-(t-p.stimDur)/p.Tri).*(1-exp(-(p.stimDur)/p.Tri));
       
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
            fi_fun2 = @(t) exp(-t/p.Tdi+500).*(fi_fun1(p.stimDur)*exp(-500) + hi_fun1(p.stimDur)/(a*p.Tdi)*(exp(a*t - 500) -exp(-500)));
        end
    end
    
    %subtract delay from ending index to not encroach on subsequent SPFR
    spfr_data.time = spfr_data.time - spfr_data.time(t_ind(1)); %because our analytical solutions assume t0 = 0

    fe_tmp = zeros(size(spfr_data.time(t_ind(1):t_ind(2)-1)));
    fe_tmp(s_ind(1):t_ind(2)-1-d) = fe_fun2(spfr_data.time(s_ind(1):t_ind(2)-1-d) - p.stimDur);
    
    fi_tmp = zeros(size(spfr_data.time(t_ind(1):t_ind(2)-1)));
    fi_tmp(t_ind(1):s_ind(1)-1) = fi_fun1(spfr_data.time(t_ind(1):s_ind(1)-1));
    fi_tmp(s_ind(1):t_ind(2)-1-d) = fi_fun2(spfr_data.time(s_ind(1):t_ind(2)-1-d) - p.stimDur);
   
    %populate a matrix with this dynamic for each stim
    fe = zeros(size(p.Ie,2),size(spfr_data.time,1));
    fi = zeros(size(p.Ii,2),size(spfr_data.time,1));
   
    %after solving for conductances at each time, we can now shift back
    %to non-delayed frame to fill vector
    t_ind = t_ind-d;
    s_ind = s_ind-d;
    for i = 1:size(Ie,1)
        fe(logical(p.Ie(i,:)),t_ind(i):(t_ind(i)+size(fe_tmp,1)-1)) = repmat(fe_tmp,1,sum(p.Ie(i,:)))';
        fi(logical(p.Ii(i,:)),t_ind(i):(t_ind(i)+size(fi_tmp,1)-1)) = repmat(fi_tmp,1,sum(p.Ii(i,:)))';
    end
    
    %caluclate conductances and steady state voltage
    ge  = p.ae*fe;
    gi  = p.ai*fi;
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = V(1:size(spfr_data.time,1));
    
else
    p.Ie = [p.Ie;zeros(2,size(p.Ie,2))];
    t_ind = [t_ind(1:end-1);s_ind(end);t_ind(end)];
    tmp = cell(length(t_ind)-1,1);
    for i = 1:length(t_ind)-1
        sol = ode45(@(t,y) vm_wrap_mov(t,y,p),[spfr_data.time(t_ind(i)) spfr_data.time(t_ind(i+1)-1)+10],y0);
        tspan = spfr_data.time(t_ind(i):t_ind(i+1)-1); % solve during bar flash
        tmp{i} = deval(sol,tspan); %solve each bit and store laterally
        gate = (p.Ie(i,:) - p.Ie(i+1,:)) > 0;
        y0 = tmp{i}(:,end); %next bit picks up where this left off
        y0(1:4:length(y0),end) = (p.b).*gate' + tmp{i}(1:4:length(y0),end).*~gate'; % set delta function, which jumps to amplitude instantly only at section that relaxes at t_ind
    end
    tmp = [tmp{:}];
    fe = tmp(2:4:length(y0),:);
    fi = tmp(4:4:length(y0),:);
    
    ge  = p.ae*fe;
    gi  = p.ai*fi;  
    
    V = ((p.Vr + ge*p.Ve + gi*p.Vi)./(1+ge+gi))';
    V = [zeros(length(spfr_data.baseSub) - length(V),1);V]; %because of delay, add bad time with zeros
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
he = y(1:4:end);
fe = y(2:4:end);
hi = y(3:4:end);
fi = y(4:4:end);

dhe = (-he') / p.Tre;
dfe = (-fe' + he')     / p.Tde;
dhi = (-hi' + p.stim*p.Ii(p.stim_num,:)) / p.Tri;
dfi = (-fi' + hi')     / p.Tdi;

tmp = [dhe;dfe;dhi;dfi];
dydt = reshape(tmp,[],1);

end
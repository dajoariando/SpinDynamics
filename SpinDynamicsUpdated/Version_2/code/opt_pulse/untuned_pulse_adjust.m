% Adjust the segment lengths of a previously-optimized pulse to cancel
% transients when fed into an untuned probe
% Written by: Soumyajit Mandal, 01/20/20
% ---------------------------------------------------------------------

function [pp_adj]=untuned_pulse_adjust(file,pulse_num)

% Load the results file
filname = file;

tmp=load(filname); results_all=tmp.results;
results=results_all{pulse_num};
siz = size(results);

if siz(2)==7 % Refocusing pulse
    params=results{5};
    sp=results{6}; pp=results{7};
    ref=1;
end

if siz(2)==10 % Excitation pulse
    params=results{8};
    sp=results{9}; pp=results{10};
    params.pcycle=0;
    ref=0;
end

sp.plt_axis=0;  sp.plt_tx=1; sp.plt_rx=0; % Set plotting parameters

% Calculate coil currrents
% ------------------------------------------------------------------------
T_90=pp.T_90; % Rectangular T_90 time
B1max=(pi/2)/(T_90*sp.gamma);
sens=sp.sens; % Coil sensitivity, T/A
amp_zero=pp.amp_zero; % Minimum amplitude for calculations

% Original pulse
% ------------------------------------------------------------------------
% Create pulse vector, including delays to allow pulse to ring down *T_90/(pi/2)
if ref % Refocusing pulse
    pp.tref=[params.tref params.tqs params.trd];
    pp.pref=[params.pref 0 0]; pp.aref=[params.aref 0 0];
    pp.Rsref=[params.Rs(2)*ones(1,length(params.tref)) params.Rs(3) params.Rs(1)];
else % Excitation pulse
    pp.tref=[params.texc params.tqs params.trd];
    pp.pref=[params.pexc 0 0]; pp.aref=[params.aexc 0 0];
    pp.Rsref=[params.Rs(2)*ones(1,length(params.texc)) params.Rs(3) params.Rs(1)];
end

[tvect,Icr,~,~] = untuned_probe_lp(sp,pp);

delt=(pi/2)*(tvect(2)-tvect(1))/T_90; % Convert to normalized time
trefc=delt*ones(1,length(tvect));
prefc=atan2(imag(Icr),real(Icr));
arefc=abs(Icr)*sens/B1max;
arefc(arefc<amp_zero)=0; % Threshold amplitude
ind=find(arefc==0); prefc(ind)=0;

if sp.plt_tx
    figure(98);
    plot(tvect/T_90,arefc.*cos(prefc)); hold on;
    plot(tvect/T_90,arefc.*sin(prefc));
    xlabel('Normalized time, t/T_{90}')
    ylabel('Normalized coil current, rotating frame')
    set(gca,'FontSize',15); set(gca,'FontWeight','bold');
end

% Modified pulse
% ------------------------------------------------------------------------
pp.N = 64; % Number of points per RF cycle
pp.NumPhases = pp.N/2; % Number of phases to quantize to

% Quantize phases
if ref
    pvec=quantize_phase(params.pref,sp,pp);
    len=length(params.pref);
    tvec=params.tref; tvec_adj=tvec;
else
    pvec=quantize_phase(params.pexc,sp,pp);
    len=length(params.pexc);
    tvec=params.texc; tvec_adj=tvec;
end

sp.plt_axis=0;  sp.plt_tx=3; sp.plt_rx=0; % Set plotting parameters

% Calculate steady-state phase shift
w=pp.w; tclk=2*pi/(w*pp.N); % Operating frequency, clock period
tau=sp.L/(sp.R+pp.Rsref(2)); % Assume Rsref is constant
theta=-atan2(w*tau,1);

% Adjust segment lengths (assume starting instant is T=0 and all sgement
% lengths are initially multiples of pi)

phi_1=pi/2-theta; del_phi=(phi_1-pvec(1)); % Amount of phase rotation
pvec=pvec+del_phi; % Rotate phases (ensures no initial transient)

for i=1:len-1 % Other segments
    alpha=mod(-(pvec(i)+pvec(i+1))/2-theta,pi);
    if alpha <= pi/2
        tadj=alpha/w;
    else
        tadj=-(pi-alpha)/w;
    end
    tadj=round(tadj/tclk)*tclk; % Quantize tadj to nearest clock edge
    tvec_adj(i)=tvec_adj(i)+tadj;
    tvec_adj(i+1)=tvec_adj(i+1)-tadj; % Adjust switching instant
end

tmp=zeros(1,2);
tmp(1)=pi/2-pvec(end)-theta; % Adjust end time to avoid final transient (try two possibilities)
tmp(2)=3*pi/2-pvec(end)-theta;
[~,ind]=min(abs(tmp)); tadj=tmp(ind)/w;
tadj=round(tadj/tclk)*tclk; % Quantize tadj to nearest clock edge
tvec_adj(end)=tvec_adj(end)+tadj;

% Create modified pulse vector, including delays to allow pulse to ring down *T_90/(pi/2)
pp_adj=pp;
pp_adj.tref=[tvec_adj params.tqs params.trd];
pp_adj.pref=[pvec 0 0];

[tvect_adj,Icr_adj,~,~] = untuned_probe_lp(sp,pp_adj);

delt=(pi/2)*(tvect_adj(2)-tvect_adj(1))/T_90; % Convert to normalized time
trefc_adj=delt*ones(1,length(tvect_adj));
prefc_adj=atan2(imag(Icr_adj),real(Icr_adj));
arefc_adj=abs(Icr_adj)*sens/B1max;
arefc_adj(arefc_adj<amp_zero)=0; % Threshold amplitude
ind=find(arefc_adj==0); prefc_adj(ind)=0;

if sp.plt_tx
    figure(99);
    plot(tvect_adj/T_90,arefc_adj.*cos(prefc_adj)); hold on;
    plot(tvect_adj/T_90,arefc_adj.*sin(prefc_adj));
    xlabel('Normalized time, t/T_{90}')
    ylabel('Normalized coil current, rotating frame')
    set(gca,'FontSize',15); set(gca,'FontWeight','bold');
end

figure(100);
if ref
    plot(params.pref+del_phi,'b--'); hold on; % Original phase
else
    plot(params.pexc+del_phi,'b--'); hold on; % Original phase
end
plot(pvec,'b-'); % Quantized phase
ylabel('Pulse phase (rad)');
xlabel('Segment number');
legend({'Original','Quantized'});
set(gca,'FontSize',15); set(gca,'FontWeight','bold');

figure(101);
plot(tvec_adj./tvec);
ylabel('Normalized segment length');
xlabel('Segment number');
set(gca,'FontSize',15); set(gca,'FontWeight','bold');
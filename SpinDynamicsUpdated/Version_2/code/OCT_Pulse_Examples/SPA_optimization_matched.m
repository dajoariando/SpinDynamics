% Summarize performance metrics of SPA-type refocusing pulses
% Assuming a matched probe

function [out] = SPA_optimization_matched

[sp,pp] = set_params_matched_SPA; % Define default system parameters

% Set excitation pulse amplitude and length
% (relative to T_90)
pp.aexc = 6; % Use aexc >> 1 to get broadband excitation
pp.texc = pp.T_90/pp.aexc; % Excitation pulse length (maintaining flip angle)
pp.tcorr=-(2/pi)*pp.texc; % Timing correction for excitation pulse

% Read in SPA refocusing pulses
[SPA_pulses]= SPA_pulse_list;
siz = size(SPA_pulses); num_pulses = siz(2);
SNR = zeros(1,num_pulses); t_p = SNR; t_E = SNR; % SPA pulses

% Define rectangular refocusing pulses
rlen = [0.6 0.8 1]; % Units of T_180
num_rpulses = length(rlen);
SNRr = zeros(1,num_rpulses); t_pr = SNRr; t_Er = SNRr;

pp.delt=0.1; % Segment length for refocusing pulsess (units of T_180)

% Rectangular
for i=1:num_rpulses
    pp.pref=[0 zeros(1,round(rlen(i)/pp.delt)) 0];
    pp.aref = [0 ones(1,length(pp.pref)-2) 0]; % Constant amplitude
    pp.tref = [pp.preDelay pp.T_180*pp.delt*ones(1,length(pp.pref)-2) pp.postDelay]; % Fixed segment length
    
    t_pr(i) = pp.delt*(length(pp.pref)-2);
    t_Er(i) = (pp.preDelay+pp.postDelay)/pp.T_180 + t_pr(i); % Calculate echo spacing (units of T_180)
    [~,~,~,~,SNRr(i)] = plot_masy_arbref_matched(sp,pp);
end

FOMtr = t_Er./(SNRr.^2); % Calculate imaging time metric (lower the better)
FOMer = t_Er.*t_pr./(SNRr.^2); % Calculate total energy metric (lower the better)

%SPA
parfor i=1:num_pulses
    pp_curr = pp;
    pp_curr.pref = [0 SPA_pulses{i} 0];
    pp_curr.aref = [0 ones(1,length(pp_curr.pref)-2) 0]; % Constant amplitude
    pp_curr.tref = [pp_curr.preDelay pp_curr.T_180*pp_curr.delt*ones(1,length(pp_curr.pref)-2) pp_curr.postDelay]; % Fixed segment length
    
    t_p(i) = pp_curr.delt*(length(pp_curr.pref)-2);
    t_E(i) = (pp_curr.preDelay+pp_curr.postDelay)/pp_curr.T_180 + t_p(i); % Calculate echo spacing (units of T_180)
    [~,~,~,~,SNR(i)] = plot_masy_arbref_matched(sp,pp_curr);
end

FOMt = t_E./(SNR.^2); % Calculate imaging time metric (lower the better)
FOMe = t_E.*t_p./(SNR.^2); % Calculate total energy metric (lower the better)

close all;

figure(1);
plot(t_pr,SNRr,'bs'); hold on;
plot(t_p,SNR,'rs-');
set(gca,'FontSize',15); set(gca,'FontWeight','bold');
xlabel('Refocusing pulse length, t/T_{180}')
ylabel('SNR')

figure(2);
plot(t_pr,FOMtr,'bo'); hold on;
plot(t_p,FOMt,'ro-');
set(gca,'FontSize',15); set(gca,'FontWeight','bold');
xlabel('Refocusing pulse length, t/T_{180}')
ylabel('FOM_{t}')

figure(3);
plot(t_pr,FOMer,'bd'); hold on;
plot(t_p,FOMe,'rd-');
set(gca,'FontSize',15); set(gca,'FontWeight','bold');
xlabel('Refocusing pulse length, t/T_{180}')
ylabel('FOM_{e}')

% Reference to T_180 rectangular pulse
t_p = [t_pr(1:end-1) t_p];
t_E = [t_Er(1:end-1) t_E];
SNR = [SNRr(1:end-1) SNR]/SNRr(end);
FOMt = [FOMtr(1:end-1) FOMt]/FOMtr(end);
FOMe = [FOMer(1:end-1) FOMe]/FOMer(end);

% Create results structure
out.t_p = t_p; out.t_E = t_E;
out.SNR = SNR; out.FOMt = FOMt; out.FOMe = FOMe;
function [SNR]=sim_inv_matched_probe_WURST_CPMG(NE,TE)

% Simulate each echo of a CPMG sequence
[sp, pp]=set_params_matched; % Define system parameters
T_90=pp.T_90; % Nominal T_90 pulse length

% Design matching network
[C1,C2]=matching_network_design2(sp.L,sp.Q,sp.f0,sp.Rs,sp.plt_mn);
sp.C1=C1; sp.C2=C2; % Save matching capacitor values

% Create (w0, w1) and coil sensitivity maps
numpts=sp.numpts; maxoffs=10;
sp.del_w=linspace(-maxoffs,maxoffs,numpts); % Linear gradient
sp.w_1=ones(1,numpts); % Uniform transmit w_1
sp.w_1r=ones(1,numpts); % Uniform receiver sensitivity

% Set acquisition parameters
tacq=(pi/2)*pp.tacq/T_90; % Normalized acquisition window length
tdw=(pi/2)*pp.tdw/T_90; % Normalized receiver dwell time
nacq=round(tacq/tdw)+1; % Number of acquired time domain points

% Set plotting parameters
sp.plt_tx = 0; 
sp.plt_rx = 0;
sp.plt_sequence = 0; % Plots /off
sp.plt_axis = 0; 
sp.plt_mn = 0; 
sp.plt_echo = 0;

% Create pulse sequence (in real time)
% --------------------------------------------
% Excitation pulse, including timing correction
[ freq,amp,phase,time,acq ] = create_WURST( 00000,2000);

texc = [time*2 2/pi 20]*pp.T_90;
pexc = [phase 0 0];
aexc = [amp 0 0];
acq_exc = [acq 0 1];
fexc = [freq 0 0];

% texc=[1 1 2/pi 1]*pp.T_90; 
% pexc=[pi/2 pi/2 0 0];
% aexc=[1 1 0 0];
% acq_exc=[0 0 0 1];
% fexc = [-10e5 10e5 0 0];

% Refocusing cycles
tref=zeros(1,3*NE); pref=tref; aref=tref; acq_ref=tref;
% for i=1:NE
%     tref((i-1)*3+1:i*3)=[(TE-pp.T_180)/2 pp.T_180 (TE-pp.T_180)/2];
%     pref((i-1)*3+1:i*3)=[0 0 0];
%     aref((i-1)*3+1:i*3)=[0 1 0];
%     acq_ref((i-1)*3+1:i*3)=[0 0 1];
% end

% Create complete pulse sequence
% pp.tp=[texc tref]; 
% pp.phi=[pexc pref];
% pp.amp=[aexc aref];
% pp.acq=[acq_exc acq_ref];

pp.tp=[texc]; 
pp.phi=[pexc];
pp.amp=[aexc];
pp.acq=[acq_exc];
pp.fp = fexc;
% Run PAP phase cycle
[~,mrx1,pnoise,f,~,~]=calc_macq_matched_probe_WURST(sp,pp);

pp.phi=[-pexc pref]; % Change phase of excitation pulse
[~,mrx2,~,~,~,~]=calc_macq_matched_probe_WURST(sp,pp);

mrx=(mrx1-mrx2); pnoise=pnoise*sqrt(2); % Apply phase cycle

figure; % Figure to show outputs
% Calculate time-domain echoes
echo_rx=zeros(NE,nacq); SNR=zeros(1,NE);
for i=1:NE
    [echo_rx(i,:),tvect]=calc_time_domain_echo_arb(mrx(i,:),sp.del_w,tacq,tdw,sp.plt_echo);
    plot((tvect*T_90/(pi/2)+i*TE)*1e3,real(echo_rx(i,:)),'b-'); hold on;
    plot((tvect*T_90/(pi/2)+i*TE)*1e3,imag(echo_rx(i,:)),'r-');
    
    SNR(i)=matched_probe_rx_snr(sp,mrx(i,:),pnoise,f);
end

set(gca,'FontSize',14);
xlabel('Time (ms)'); ylabel('Received echoes');

end


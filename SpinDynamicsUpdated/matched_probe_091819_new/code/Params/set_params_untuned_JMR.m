%SET_PARAMS_UNTUNED Summary of this function goes here
%   Detailed explanation goes here

function [sp,pp] = set_params_untuned_JMR
% Global parameters
% --------------------------------------------
sp.k = 1.381e-23; % J/K
sp.T = 300; % Sample temperature

% System parameters
% --------------------------------------------
sp.f0 = 0.5e6; % Target matching frequency (= Larmor frequency), Hz
sp.fin = 0.5e6; % Input frequency, Hz
sp.w0 = 2*pi*sp.fin;
pp.w = sp.w0;

% Coil parameters
% --------------------------------------------
sp.L = 10e-6; % H
sp.Q = 50;
sp.R = 2*pi*sp.f0*sp.L/sp.Q;

% Tuning Params
%-------------------------------------
sp.C = 1/((2*pi*10*sp.f0)^2*sp.L); % Parasitic cap, fsr = 10 x f0

% Transmitter parameters
% --------------------------------------------
sp.Rs = 2; % Series resistance, Ohms

% Receiver parameters
% --------------------------------------------
sp.Rin = 1e6; % Input impedance, Ohms
sp.Cin = 5e-12;
sp.Rd = 1e6; % Damping resistance
sp.Rdup = 0.2; % Duplexer on resistance

% Input transformer parameters (Coilcraft PWB series)
sp.Nrx = 4; sp.krx = 0.9996; % Turns ratio, coupling factor
sp.L1 = 75e-6; sp.R1 = 0.26; % Transformer primary
sp.L2 = 1250e-6; sp.R2 = 0.91; % Transformer secondary

sp.NF = 1; % Noise figure, dB
sp.vn = 0.45e-9; 
sp.in = 2e-15; % Receiver input noise voltage [V/sqrt(Hz)] and current [A/sqrt(Hz)]

% Simulation parameters
% --------------------------------------------
sp.m0=1; % Initial magnetization
sp.mth=1; % Asymptotic / thermal magnetization
sp.numpts=2000;  
sp.maxoffs=10;
sp.del_w=linspace(-sp.maxoffs,sp.maxoffs,sp.numpts); % Static Gradient

% Matched filter type
sp.mf_type=2; % 1 -> matched (white noise), 2 -> matched (colored noise)

% Plotting parameters
% --------------------------------------------
sp.plt_tx = 0; 
sp.plt_rx = 1; 
sp.plt_sequence = 1; % Plot on/off
sp.plt_axis = 1; 
sp.plt_mn = 1; 
sp.plt_echo = 1;

% Pulse sequence parameters
% --------------------------------------------
pp.N=32; % quantization step = N x input RF frequency
pp.T_90=25e-6;
pp.T_180=2*pp.T_90; % Rectangular T_90 and T_180
pp.psi=0; % Absolute RF phase at t=0
pp.preDelay = 20e-6;
pp.postDelay = 50e-6;

% Excitation pulse
pp.texc=[1]*pp.T_90;
pp.pexc=[pi/2]; 
pp.aexc=[1];
pp.tcorr=-(2/pi)*pp.T_90; % Timing correction for excitation pulse

% Refocusing cycle
pp.tref=[pp.preDelay pp.T_180 pp.postDelay]; 
pp.pref=[0 0 0]; 
pp.aref=[0 1 0];
pp.Rsref=[2 2 20];

pp.tacq=[5]*pp.T_180; % Acquisition time for observing echo
pp.tdw=0.5e-6; % Receiver dwell time

pp.amp_zero=1e-4; % Minimum amplitude for calculations
end
% Use Martin/Yi-Qiao's terminology to simulate dynamics of spin-1/2
% ensembles. Pulses have arbitrary power levels.
% ------------------------------------
% params structure:
% ------------------------------------
% tp = RF pulse durations in sec (normalized  to w1 = 1)
% phi = RF pulse phases in radians
% amp = pulse amplitude (zero for free precession)
% grad = gradient strength (zero for no gradient)
% acq = acquire signal if 1
% len_acq = duration for observing echo in sec (normalized  to w1 = 1)
% del_w, w_1 = (del_w0, w1) map
% del_wg = (change in del_w0) map due to a gradient pulse (grad = 1)
% T1n, T2n = relaxation time constant map (normalized  to w1 = 1)
% All coherence pathways are considered in this simulation
% ------------------------------------
% Soumyajit Mandal 08/26/10
% Allowed arbitrary pulse amplitudes 02/25/11
% Use normalized w1, allow relaxation during free precession 03/16/11
% Use 2rd order Trotter approximatation for relaxation during pulses 10/01/11
% Use Rodrigues' rotation formula to speed up simulation (~4x) 10/04/11
% Allow arbitrary w0, w1 maps 10/04/11
% Use completely vectorized code for further speed (~3x speedup): 05/06/14
% Reintroduce relaxation during free precession: 09/09/19
% Add separate input variable del_wg for changes in del_w (useful for imaging):
% 09/09/19

function [macq]=sim_spin_dynamics_arb7(params)

% Load parameters
tp=params.tp;
phi=params.phi;
amp=params.amp; 
acq=params.acq;
grad=params.grad;
len_acq=params.len_acq; 
del_w=params.del_w;
del_wg=params.del_w;
w_1=params.w_1;
T1n=params.T1n;
T2n=params.T2n;
mth = params.mth;

% Normalize tp


m0=params.m0; % Initial magnetization vector amplitude
% mth=sp.mth; % Thermal magnetization vector amplitude

numpts=length(del_w);
if length(w_1)~=numpts
    disp('Error: w0 and w1 vectors have different lengths!');
    return;
end
if length(del_wg)~=numpts
    disp('Error: w0 and w0g vectors have different lengths!');
    return;
end
% window = sinc(del_w*len_acq/(2*pi)); % window function for acquisition
window = sinc(del_w*1/(2*pi)); % window function for acquisition
window = window./sum(window);

% Calculate magnetization spectrum (no diffusion)
mvect=zeros(3,numpts); % Magnetization vectors
mvect(1,:)=m0; % Initial vectors are along z-axis

acq_cnt=0; % Acquisition counter
nacq=sum(acq); % Number of acquisitions
macq=zeros(nacq,numpts);

% Evolution of magnetization
num_pulses=length(phi); del_w0=del_w;
for j=1:num_pulses
    del_w=del_w0+grad(j)*del_wg; % Apply gradient pulse
    
    if amp(j)>0
        w1=amp(j)*w_1;
        Omega=sqrt(w1.*w1+del_w.*del_w);
        mat=calc_matrix_elements(del_w,w1,Omega,tp(j),phi(j)); % RF pulses
        mlong=zeros(1,numpts); % No longitudinal relaxtion
    else
        mat=calc_fp_matrix_elements(del_w,tp(j),T1n,T2n); % Free precession
        mlong=mth.*(1-exp(-tp(j)./T1n)); % Longitudinal relaxation
    end
    
    tmp=mvect;
    mvect(1,:)=mat.R_00.*tmp(1,:)+mat.R_0m.*tmp(2,:)+mat.R_0p.*tmp(3,:)+mlong; % M0
    mvect(2,:)=mat.R_m0.*tmp(1,:)+mat.R_mm.*tmp(2,:)+mat.R_mp.*tmp(3,:); % M-
    mvect(3,:)=mat.R_p0.*tmp(1,:)+mat.R_pm.*tmp(2,:)+mat.R_pp.*tmp(3,:); % M+
    
    if acq(j) % Acquire spectrum (-1 coherence)
        acq_cnt=acq_cnt+1;
        macq(acq_cnt,:)=conv(mvect(2,:),window,'same'); % Convolve with acquisition window
    end
end

% Calculate matrix elements for RF pulses, neglect relaxation
function R = calc_matrix_elements(del_w,w1,Omega,tp,phi)

dw=del_w./Omega; 
dw_2=dw.*dw; 
w1n=w1./Omega; 
w1n_2=w1n.*w1n;
ph=exp(1i*phi); 
sn=sin(Omega*tp);
cs=cos(Omega*tp);

R.R_00=dw_2+w1n_2.*cs;
R.R_0p=0.5*w1n.*(dw.*(1-cs)-1i*sn)*conj(ph); R.R_0m=conj(R.R_0p);
R.R_p0=w1n.*(dw.*(1-cs)-1i*sn)*ph; R.R_m0=conj(R.R_p0);
R.R_pp=0.5*(w1n_2+(1+dw_2).*cs)+1i*dw.*sn; R.R_mm=conj(R.R_pp);
R.R_pm=0.5*w1n_2.*(1-cs)*ph*ph; R.R_mp=conj(R.R_pm);

% For free precession, w1 = 0, include relaxation
function R = calc_fp_matrix_elements(del_w,tf,T1n,T2n)

numpts=length(del_w);
R.R_00=exp(-tf./T1n).*ones(1,numpts);
R.R_0p=zeros(1,numpts); R.R_0m=conj(R.R_0p);
R.R_p0=zeros(1,numpts); R.R_m0=conj(R.R_p0);
R.R_pp=exp(-tf./T2n).*exp(1i*del_w*tf); R.R_mm=conj(R.R_pp);
R.R_pm=zeros(1,numpts); R.R_mp=conj(R.R_pm);
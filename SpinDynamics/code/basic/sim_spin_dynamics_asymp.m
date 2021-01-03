% Use Martin/Yi-Qiao's terminology to simulate dynamics of spin-1/2
% ensembles. Find asymptotic echo shape from effective rotation axis.
% ------------------------------------
% T_90 = duration of 90 degree RF pulse
% tp = RF pulse durations in us
% phi = RF pulse phases in radians
% tf = spacing between pulses in us
% neff = effective rotation axis
% all coherence pathways are considered in this simulation
% all RF pulses are assumed to have the same power level
% ------------------------------------
% Soumyajit Mandal 08/26/10

function [echo,tvect]=sim_spin_dynamics_asymp(T_90,T_180,tp,phi,tf,neff)

% Units conversions
T_90=T_90/1e6;
T_180=T_180/1e6;
tp=tp/1e6;
tf=tf/1e6;
i=sqrt(-1);

m0=1; % Initial magnetization vector amplitude

% Uniform distribution of del_w0 (uniform gradient)
maxoffs=10;
numpts=1e3; % number of spin vectors to simulate
necho=6e2; % number of time points for observing echo
dt=1e-6; % sampling period for observing echo

% Check if neff has the right length
siz_neff=size(neff);
if siz_neff(3)~=numpts
    disp('Error! Effective rotation axis vector of the wrong length!');
    echo=0; tvect=0;
    return;
end

delt=dt*necho; % Final time period for observing echo = 600 us (default)
tvect=linspace(-delt/2,delt/2,necho);

w1=pi/(2*T_90);
del_w=(pi/T_180)*linspace(-maxoffs,maxoffs,numpts);

Omega=sqrt(w1^2+del_w.^2);

% Calculate echo spectrum (no diffusion)
mvect=zeros(3,1,numpts); % Magnetization vectors
mvect(1,1,:)=m0*ones(1,1,numpts); % Initial mag vectors are along z-axis

num_pulses=length(phi);
for j=1:num_pulses
    mat=calc_matrix_elements(del_w,w1,Omega,tp(j),phi(j)); % RF pulses
    for k=1:numpts
        mvect(:,:,k)=mat(:,:,k)*mvect(:,:,k);
    end
    
    if tf(j)>0
        mat=calc_fp_matrix_elements(del_w,tf(j)); % Free precession
        for k=1:numpts
            mvect(:,:,k)=mat(:,:,k)*mvect(:,:,k);
        end
    end
end

% Asymptotic magnetization is spin-locked to the effective rotation axis
for k=1:numpts
   mvect(:,1,k)=dot(neff(:,1,k),mvect(:,1,k))*neff(:,1,k); 
end

% Only -1 coherence is visible
trans=zeros(1,numpts);

for k=1:numpts
    trans(k)=mvect(map_pi(-1),1,k);
end

% Calculate time-domain waveform
echo=zeros(1,necho);

for j=1:necho
    echo(j)=sum(trans.*exp(-i*del_w*tvect(j)))/numpts; % -1 coherence at time of detection
end

% Mapping from coherence to index
% [0,-1,+1] -> [1,2,3]
function [ind]=map_pi(coh)

if coh==0
    ind=1;
else if coh==-1
        ind=2;
    else
        ind=3;
    end
end

% Calculate matrix elements for RF pulses, neglect relaxation
function R = calc_matrix_elements(del_w,w1,Omega,tp,phi)

i=sqrt(-1);
R_pp=0.5*((w1./Omega).^2+(1+(del_w./Omega).^2).*cos(Omega*tp))+i*(del_w./Omega).*sin(Omega*tp);
R_mm=conj(R_pp);
R_00=(del_w./Omega).^2+(w1./Omega).^2.*cos(Omega*tp);
R_p0=(w1./Omega).*((del_w./Omega).*(1-cos(Omega*tp))-i*sin(Omega*tp))*exp(i*phi);
R_m0=conj(R_p0);
R_0p=0.5*(w1./Omega).*((del_w./Omega).*(1-cos(Omega*tp))-i*sin(Omega*tp))*exp(-i*phi);
R_0m=conj(R_0p);
R_pm=0.5*(w1./Omega).^2.*(1-cos(Omega*tp))*exp(i*2*phi);
R_mp=conj(R_pm);

% Rotation matrix
R=zeros(3,3,length(del_w));
R(1,1,:)=R_00; R(1,2,:)=R_0m; R(1,3,:)=R_0p;
R(2,1,:)=R_m0; R(2,2,:)=R_mm; R(2,3,:)=R_mp;
R(3,1,:)=R_p0; R(3,2,:)=R_pm; R(3,3,:)=R_pp;

% For free precession, w1 = 0, neglecting relaxation
function R = calc_fp_matrix_elements(del_w,tf)

i=sqrt(-1);
R_pp=cos(del_w*tf)+i*sin(del_w*tf);
R_mm=conj(R_pp);
R_00=ones(1,length(del_w));
R_p0=zeros(1,length(del_w));
R_m0=conj(R_p0);
R_0p=zeros(1,length(del_w));
R_0m=conj(R_0p);
R_pm=zeros(1,length(del_w));
R_mp=conj(R_pm);

% Rotation matrix
R=zeros(3,3,length(del_w));
R(1,1,:)=R_00; R(1,2,:)=R_0m; R(1,3,:)=R_0p;
R(2,1,:)=R_m0; R(2,2,:)=R_mm; R(2,3,:)=R_mp;
R(3,1,:)=R_p0; R(3,2,:)=R_pm; R(3,3,:)=R_pp;
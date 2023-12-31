% Calculate receiver filtering function for an untuned probe (assuming an
% input transformer is present)

function [mrx,SNR,tf] = untuned_probe_rx(sp,pp,macq)

k=sp.k; T=sp.T;
L=sp.L; R=sp.R; C=sp.C;
Cin=sp.Cin; Rin=sp.Rin; Rd=sp.Rd;
Rdup=sp.Rdup;

Nrx=sp.Nrx; krx=sp.krx;
L1=sp.L1; R1=sp.R1; % Transformer primary
L2=L1*Nrx^2; R2=Nrx*R1; % Transformer secondary

vn=sp.vn; in=sp.in;

mf_type=sp.mf_type;

w0=sp.w0; % Larmor frequency
del_w=sp.del_w;
w1_max=(pi/2)/pp.T_90;

s=1i*(w0+del_w*w1_max); f=imag(s)/(2*pi); % Un-normalized frequency axis
Yin=s*Cin+1/Rin+1/Rd; % Input admittance
Yp=(s*C+1./(s*L+R)); Zp=1./Yp; % Impedance of parallel resonator

Nv=krx*Nrx*L1/(L+L1); % Voltage gain ratio
Zeff=(Zp+Rdup+R1)*Nv^2+s*L2*(1-krx^2)+R2; % Effective source impedance

tf=Nv./(1+Zeff.*Yin);
Zs=Zeff./(1+Yin.*Zeff);

% Calculate receiver waveform and noise PSD
mrx=macq.*tf.*(2*pi*f/w0).^2;
% Note multiplication of mrx(w) by w^2 to account for polarization and
% inductive detection efficiency factors

vni2=4*k*T*R*abs(tf).^2;
pnoise=vn.^2+in.^2*abs(Zs).^2+4*k*T*(Rdup+R1+R2/Nv^2)*abs(tf).^2+vni2;

% Matched filtering (maximize SNR at t=0)
switch mf_type
    case 0 % Rectangular in time over window
        theta=atan2(sum(imag(mrx)),sum(real(mrx)));
        %T_W=pp.len_acq; % window length = len_acq
        T_W=0.8*pi; % window length = 0.8*T_180
        mf=sinc(del_w*T_W/(2*pi))*exp(-1i*theta); 
    case 1 % Matched filter (for white noise)
        mf=conj(mrx);
    case 2 % Matched filter (for colored noise)
        mf=conj(mrx)./pnoise;
end
mf=mf/sqrt(trapz(del_w,abs(mf).^2)); % Normalize matched filter amplitude

vsig=trapz(del_w,mrx.*mf); % Signal amplitude (at t=0)
vnoise=sqrt(trapz(f,pnoise.*abs(mf).^2)); % Total integrated noise (V_rms)
SNR=real(vsig)/vnoise; % Signal-to-noise ratio (SNR) per scan in voltage units

NF=pnoise./vni2; % Calculate receiver noise figure

if sp.plt_rx
    figure(6);
    plot(del_w,real(macq),'b--','LineWidth',1); hold on;
    plot(del_w,imag(macq),'r--','LineWidth',1);
    plot(del_w,real(mrx)/max(abs(tf)),'b-','LineWidth',1); hold on;
    plot(del_w,imag(mrx)/max(abs(tf)),'r-','LineWidth',1);
    ylabel('S_{asy}(\omega), S_{rx}(\omega)');
%    ylabel('S_{rx}(\omega)');
    xlabel('\Delta\omega_{0}/\omega_{1,max}');
    set(gca,'FontSize',15); set(gca,'FontWeight','bold');
    
    figure(7);
    subplot(1,2,1);
    plot(del_w,pnoise/vn^2,'k-'); hold on;
    ylabel('N(\omega)');
    xlabel('\Delta\omega_{0}/\omega_{1,max}');
    title('Filter functions');
    set(gca,'FontSize',15); set(gca,'FontWeight','bold');
    
    subplot(1,2,2);
    plot(del_w,abs(mf),'r-'); hold on;
    ylabel('|H(\omega)|');
    xlabel('\Delta\omega_{0}/\omega_{1,max}');
    set(gca,'FontSize',15); set(gca,'FontWeight','bold');
     
    figure(8);
    subplot(2,1,1); plot(del_w,abs(tf),'b-');
    ylabel('mag(TF)');
    title('Receiver gain A(\omega)');
    set(gca,'FontSize',15); set(gca,'FontWeight','bold');
    
    subplot(2,1,2); plot(del_w,(180/pi)*angle(tf),'r-');
    ylabel('phase(TF)');
    xlabel('\Delta\omega_{0}/\omega_{1,max}');
    set(gca,'FontSize',15); set(gca,'FontWeight','bold');
    
    figure(9);
    plot(del_w,10*log10(NF),'b-'); hold on;
    ylabel('NF (dB)');
    xlabel('\Delta\omega_{0}/\omega_{1,max}');
    title('Receiver NF');
    set(gca,'FontSize',15); set(gca,'FontWeight','bold');
end

% mrx=mrx/max(abs(tf)); % Normalize amplitude of received magnetization
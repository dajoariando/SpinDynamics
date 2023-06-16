% Create typical w0 and w1 fields for a single-sided sensor
% Static gradient: y (assumed uniform) (T/m)
% FOV: (x,z) plane (m), uniform field

function [sp]=create_fields_single_sided(sp)
% settings
load_ext_graddata = 1; % LOAD EXTERNAL GRADIENT PLANE DATA
force_flat_tx_plane = 0; % FORCE THE TX PLANE TO BE FLAT, IGNORE load_ext_txdata below
load_ext_txdata = 1; % LOAD EXTERNAL TRANSMIT PLANE DATA


% Fixed offset: strong gradient in y-direction
ny=sp.ny; maxoffs=sp.maxoffs; % y-axis: ny is npoints in y direction. maxoffs is max offset in magnetic field (due to parasitic gradient)
del_w0y=linspace(-maxoffs,maxoffs,ny); % Linear y-gradient

if sp.plt_fields
    figure;
    plot(del_w0y); title('\Delta\omega_{0} along y')
end

% Pulsed offset: gradients in (x,z) plane
nx=sp.nx; nz=sp.nz; % Size of (x,z) plane
x1=linspace(1,nx,nx); z1=linspace(1,nz,nz);
del_wgx=linspace(-1,1,nx); % Normalized to +/- 1
del_wgz=linspace(-1,1,nz);

% GRADIENT PLANE DEFINITION
if (load_ext_graddata)
    % load external gradient field (david ariando 230616)
    gx_file = 'coil_gx_plane.txt'; % load files
    gz_file = 'coil_gz_plane.txt'; % load files
    [~,~,del_wgx] = get_interp_field_plane_data  (gx_file, nx, nz);
    [~,~,del_wgz] = get_interp_field_plane_data  (gz_file, nx, nz);
    % normalize to (-1,1)
    del_wgx=del_wgx./max(del_wgx,[],"all");
    del_wgz=del_wgz./max(del_wgz,[],"all");
    % format the data to the original data provided by Dr. Mandal
    del_wgx = del_wgx';
    del_wgz = del_wgz' .* -1;
    % Reshape to vectors
    del_wgxv=reshape(del_wgx,1,nx*nz);
    del_wgzv=reshape(del_wgz,1,nx*nz);
else
    % set the 2D-space gradient (original by Dr. Mandal)
    del_wgx=repmat(del_wgx',1,nz); %%% create nz*nx matrix of gradient to represent 2D-space gradient
    del_wgz=repmat(del_wgz,nx,1); %%% create nz*nx matrix of gradient to represent 2D-space gradient
    % Reshape to vectors %%% possibly for repetitive calculation?
    del_wgxv=reshape(del_wgx,1,nx*nz);
    del_wgzv=reshape(del_wgz,1,nx*nz);
end


if sp.plt_fields
    figure;
    subplot(1,2,1); imagesc(z1,x1,del_wgx); colorbar; title('G_x');
    subplot(1,2,2); imagesc(z1,x1,del_wgz); colorbar; title('G_z');
end

% TRANSMIT COIL PLANE DEFINITION
if (load_ext_txdata)
    % load external gradient field (david ariando 230616)
    tx_file = 'coil_tx_plane.txt'; % load files
    [~,~,xz] = get_interp_field_plane_data  (tx_file, nx, nz);
    % normalize to (-1,1)
    xz=xz./mean(xz,"all");
    % format the data to the original data provided by Dr. Mandal
    xz = xz';
    % Reshape to vectors
    xzv=reshape(xz,1,nx*nz);
    xzv=reshape(xz,1,nx*nz);
else % the original gaussian RF transmit
    % Symmetric Gaussian w1 in (x,z) plane, uniform in y
    sigma_x=1*nx; sigma_z=1*nz; % Width of gaussian (adjust as needed)
    [Z1,X1]=meshgrid(z1,x1); XZ=[X1(:) Z1(:)];
    xz=mvnpdf(XZ,[round(nx/2),round(nz/2)],[sigma_x^2,0;0,sigma_z^2]);
    xz=reshape(xz,nx,nz);
    xz=xz/max(max(xz)); % Normalize to value at center
    xzv=reshape(xz,1,nx*nz); % Reshape to vector
end
if force_flat_tx_plane % force tx coil plane to be flat 1
    xz=ones(nx,nz);
    xzv=ones(1,nx*nz);
end

if sp.plt_fields
    figure;
    if nx>1 && nz>1
        surf(z1,x1,xz); title('Normalized RF amplitude in (x,z) plane');
    end
end

% Sample properties
rho=sp.rho; T1map=sp.T1map; T2map=sp.T2map;
rhov=reshape(rho,1,nx*nz); % Reshape to vector
T1mapv=reshape(T1map,1,nx*nz);
T2mapv=reshape(T2map,1,nx*nz);

if sp.plt_fields
    figure;
    subplot(1,3,1); imagesc(z1,x1,rho); title('Spin density')
    subplot(1,3,2); imagesc(z1,x1,T1map); title('T1')
    subplot(1,3,3); imagesc(z1,x1,T2map); title('T2')
end

% Create offset frequency, RF amplitude, and sample vectors
% Created by flattening (x,z) for each y to a 1D list, then repeating for each y
sp.numpts=ny*nx*nz;
del_w0=zeros(1,sp.numpts); w_1=del_w0; del_wx=del_w0; del_wz=del_w0;
m0=del_w0; mth=del_w0; T1=del_w0; T2=del_w0;
for i=1:sp.ny
    del_w0((i-1)*nx*nz+1:i*nx*nz)=del_w0y(i)*ones(1,nx*nz);
    del_wx((i-1)*nx*nz+1:i*nx*nz)=del_wgxv;
    del_wz((i-1)*nx*nz+1:i*nx*nz)=del_wgzv;
    w_1((i-1)*nx*nz+1:i*nx*nz)=xzv;
    m0((i-1)*nx*nz+1:i*nx*nz)=rhov;
    mth((i-1)*nx*nz+1:i*nx*nz)=rhov;
    T1((i-1)*nx*nz+1:i*nx*nz)=T1mapv;
    T2((i-1)*nx*nz+1:i*nx*nz)=T2mapv;
end

sp.del_w=del_w0;
sp.w_1=w_1;
sp.del_wx=del_wx;
sp.del_wz=del_wz;
sp.m0=m0; sp.mth=mth; sp.T1=T1; sp.T2=T2;

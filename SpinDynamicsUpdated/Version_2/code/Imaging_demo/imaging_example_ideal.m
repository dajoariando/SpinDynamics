% Imaging example
% ----------------------------------------------
%parpool('local',32)
% ----------------------------------------------
% Define parameters
% ----------------------------------------------

% add necessary path
addpath(genpath('../.'));

% Pulse sequence
% ----------------------------------------------
params.NE=6; % Number of echoes
params.TE=0.2e-3; % Echo period (sec)
params.Tgrad=0.5e-3; % Gradient length (sec)
pixels = 64;
FOV = 64;
load_ext_image = 0;

if load_ext_image
    % load external image
    flower = imread('../Images/flower.png');
    flowerResize = imresize(flower,[pixels pixels]);
    IM = rgb2gray(flowerResize);
else % create circle sensitivity map
    % use a simple circle map
    IM = create_circle (pixels, pixels, pixels/2*0.95);
    
    % or use concentric circle map
    a = create_circle (pixels, pixels, pixels/2*0.95);
    b = create_circle (pixels, pixels, pixels/2*0.85);
    c = create_circle (pixels, pixels, pixels/2*0.75);
    d = create_circle (pixels, pixels, pixels/2*0.65);
    e = create_circle (pixels, pixels, pixels/2*0.55);
    f = create_circle (pixels, pixels, pixels/2*0.45);
    g = create_circle (pixels, pixels, pixels/2*0.35);
    h = create_circle (pixels, pixels, pixels/2*0.25);
    i = create_circle (pixels, pixels, pixels/2*0.15);
    IM = a-b+c-d+e-f+g-h+i;
end




% Sample parameters: change as needed to get interesting images
% ----------------------------------------------
%params.rho=ones(16,16); % Spin density map (kind of boring right now)
params.rho = IM;
params.T1map=5e-3*ones(pixels,pixels); % T1 map (also boring)
params.T2map=5e-3*ones(pixels,pixels); % T2 map (also boring)

%params.T2map = IM;
% Image parameters
% ----------------------------------------------
params.pxz=[pixels,pixels]; % Image size in pixels (x,z)
params.FOV=[FOV,FOV]; % FOV in pixel units (x,z)

% ----------------------------------------------
% Run simulation
% ----------------------------------------------
[echo_int_all]=sim_cpmg_ideal_probe_img(params);

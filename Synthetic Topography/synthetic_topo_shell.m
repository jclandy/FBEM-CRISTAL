function [PosT_s,PosT_si,surface_type] = synthetic_topo_shell(op_mode,topo_type,pitch,roll,sigma_surf_s,sigma_surf_si,l_surf,H_surf,dx,L_w,L_h,D_off,f_p)

%% Generates xyz vectors of synthetic sea ice topography

% Virtual sea ice surface topography with Gaussian, lognormal or fractal
% roughness properties

% Input:
% op_mode: 1 = pulse-limited, 2 = SAR
% topo_type: 1 = Gaussian, 2 = lognormal, 3 = fractal
% sigma_surf = rms roughness height
% l_surf = correlation length
% H_surf = Hurst parameter for fractal surfaces
% dx = grid resolution
% L_w = lead width
% L_h = lead depth
% D_off = lead distance off nadir
% f_p = melt pond fraction

% Output:
% PosT = n x 3 matrix of the xyz surface topography coordinate vertices
% surface_type: 0 = lead/ocean, 1 = sea ice

% (C) Jack Landy, University of Bristol, 2018

warning('off','all')

%% Create grid

% (Can be modified but these are default grid sizes)
if op_mode == 1
    L = 8000; % across-track diameter of grid, m
    W = 8000; % along-track diameter of grid, m
else
    L = 8000; % across-track diameter of grid, m
    W = 400; % along-track diameter of grid, m
end

[x,y] = meshgrid(-W/2+dx/2:dx:W/2-dx/2,-L/2+dx/2:dx:L/2-dx/2);

%% Generate topography

if topo_type == 1 % Gaussian
    [z_s,z_si,~,~] = rsgene2D_anisotrop(W/dx,L/dx,W,L,sigma_surf_s,sigma_surf_si,l_surf,l_surf,0);
elseif topo_type == 2 % Lognormal
    [z_s,z_si,~,~] = rsgene2D_anisotrop(W/dx,L/dx,W,L,sigma_surf_s,sigma_surf_si,l_surf,l_surf,1);
else % Fractal
    % [z_s,z_si,~,~] = artificial_surf(sigma_surf_si,H_surf,W,W/dx,L/dx,(2*pi)/(l_surf*dx));
end

%% Add lead

surface_type = ones(size(z_si)); % sea ice = 1, lead/ocean = 0

if L_w>0
    z_si(L/(dx*2)-0.5*L_w/dx-D_off/dx:L/(dx*2)+0.5*L_w/dx-D_off/dx,:) = -L_h;
    surface_type(L/(dx*2)-0.5*L_w/dx-D_off/dx:L/(dx*2)+0.5*L_w/dx-D_off/dx,:) = 0;
else
end

%% Add melt ponds

if f_p>0
    [z_si,surface_type] = add_melt_ponds(z_si,surface_type,f_p); % topo still referenced to mean height
else
end

%% Finalize

PosT_s = [x(:) y(:) z_s(:)];

PosT_si = [x(:) y(:) z_si(:)];
surface_type = surface_type(:);

warning('on','all')

end


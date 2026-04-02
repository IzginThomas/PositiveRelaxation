Nx = 100;
c = 1;



% Links- und rechtsseitige Dichte (rhol,rhor)

rhol = 0.8;
rhor = 1;

% Links- und rechtsseitiger Geschwindigkeit (vl, vr)

rhovl =  1e-3;
rhovr =  1e-2;

xstart = 0;
xend = 1;

xspan = xend - xstart;
bc = 'periodic';

x = linspace(xstart,xend,Nx);
dx = xspan/ (Nx - 1);

tstart = 0;
tend = 1;
ylabels = {'$\rho$', '$\rho v$'};



% Membranposition

membran_position = 0.5;

l_num_of_cells = round(Nx * (membran_position - xstart) / xspan);

ind_rhol = 1:l_num_of_cells;
ind_rhor = l_num_of_cells + 1: Nx;
ind_rhovl = Nx + 1 : Nx + l_num_of_cells;
ind_rhovr = Nx + l_num_of_cells + 1: 2*Nx;

y0(ind_rhol) = rhol;
y0(ind_rhor) = rhor;
y0(ind_rhovl) = rhovl;
y0(ind_rhovr) = rhovr;
y0 = y0';

pos_ind = 1:Nx;
inds = [pos_ind; setdiff(1:2*Nx,pos_ind)];
% exp_ind = setdiff(1:2*Nx,pos_ind);

% (Initial) Time step size
para = c;
CFL = 1;
dt0 = CFL * dx/(max(c + abs(y0(Nx+1:end)))); % CFL
dt0 =  1e-2;


E = [];

phys_flux = @(y,para) [y(Nx+1:2*Nx); y(Nx+1:2*Nx)./y(1:Nx) + para^2*y(1:Nx)];

num_flux = @euler_flux;

PD = @PD_isothEuler;

pde_plot = true;

% Relaxation
entropy_flag = true;
if entropy_flag
    % Entropy
    eta = str2func(['@(Y,para) sum( (Y(end/2 + 1 : end, :)).^2 ./ (2*Y(1:end/2,:)) + Y(1:end/2,:) .* para^2 .* log(Y(1:end/2,:)) )*' num2str(dx)]);
    entropy_cons = true;
    % Entropy prime
    eta_prime =str2func(['@(y,para) [-(y(end/2 + 1 : end)./y(1:end/2)).^2 + para^2.*(log(y(1:end/2)) + 1); y(end/2 + 1 : end)./y(1:end/2)] *' num2str(dx)]);
end
% Linear advection:  u_t + para*u_x=0

% Right hand side
PD = @PD_scalar_pde;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time interval
tstart = 0;
tend = 1;

% Spatial domain
Nx = 100;
xstart = 0;
xend = 2;
xspan = xend - xstart;
bc = 'periodic';

x = linspace(xstart,xend,Nx);
% x(end)=[];
% Nx=Nx-1;
dx = x(2)-x(1);


% (Initial) Time step size
para = 1;
CFL = 1;
dt0 = CFL * dx/para; % CFL

% Initial values
%y0_func =@(x) 1.*(x<=0.5) + (1.5-x + 1e-2).*((x>0.5) & (x<1.5)) + 1e-2.*(x>=1.5) ; %sin(2*pi*x) + 1;
y0_func = @(x) 1.9*sin(pi*x) + 2;
y0 = (y0_func(x))';
 % exsol = @(t,x) y0_func(mod(x-t*para,2));

% Compound composition matrix
E=[];

% fluxes

phys_flux = @(y,para) para*y;

num_flux = @harm_mean;

pde_plot = true;

% Relaxation
entropy_flag = true;
if entropy_flag
    % Entropy
    switch func2str(num_flux)
        case  'log_mean'
            eta = str2func(['@(Y) sum(Y .* log(Y) - Y)*' num2str(dx)]);
            % Entropy prime
            eta_prime =str2func(['@(y) (log(y)) *' num2str(dx)]);
        case 'geo_mean'
            eta = str2func(['@(Y) -sum(sqrt(Y))*' num2str(dx)]);
            % Entropy prime
            eta_prime =str2func(['@(y) -(1./sqrt(y)) *' num2str(dx)]);
        case 'harm_mean'
            eta = str2func(['@(Y) -sum(1./Y)*' num2str(dx)]);
            % Entropy prime
            eta_prime =str2func(['@(y) (1./(y.^2)) *' num2str(dx)]);
        case 'arith_mean'
            eta = str2func(['@(Y) sum(Y.^2/2)*' num2str(dx)]);
            % Entropy prime
            eta_prime =str2func(['@(y) y *' num2str(dx)]);
            %mass_cons_relax = false; 
    end
    entropy_cons = true;

end
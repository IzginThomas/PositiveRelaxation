% Right hand side
PD = @PD_PME;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time interval
tstart = 1;
tend = 2;

% Spatial domain
Nx = 160;
xstart = -6;
xend = 6;

x = linspace(xstart,xend,Nx);
dx = (xend - xstart) / (Nx - 1);

% (Initial) Time step size
dt0 = dx;
% dt0 = 20*dx;

% Initial values
para = 5;


%exact solution
exsol = @(t,x,m) t.^(-1/(m+1)).*max(0, 1- (m-1).*abs(x).^2./((m+1).*2.*m.*t.^(2./(m+1)))).^(1./(m-1));
m = para;
exsol = @(t,x) exsol(t,x,m);

y0 = (exsol(1,x))';

% Compound composition matrix
E=[];

pde_plot = true;

% Relaxation
entropy_flag = true;
if entropy_flag
    % Entropy
    % Entropy
    eta = str2func(['@(Y) sum(Y.^2)/2*' num2str(dx^2)]);
    entropy_cons = false;
    % Entropy prime
    eta_prime =str2func(['@(y) y*' num2str(dx^2)]);
end
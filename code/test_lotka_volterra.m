%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Right hand side
PD = @PD_lotka_volterra;

% Time interval
tstart = 0; 
tend = 100;

% Initial value

y0 = [2;2];

% Time step size
dt0 = 1;

E = [];

% Relaxation
entropy_flag = true;
if entropy_flag
        % Entropy
    eta = @(Y) -( log(Y(1,:))-Y(1,:) + 2*log(Y(2,:)) - Y(2,:) );
    % para = scaleode;

    entropy_cons = true;
    % Entropy prime
    eta_prime = @(y) [-1/y(1) + 1; -2/y(2) + 1];
    mass_cons_relax = false;
end
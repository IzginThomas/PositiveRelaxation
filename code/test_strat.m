%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stratospheric reaction problem                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Right hand side
PD = @PD_strat;

% Time interval
tstart = 12;%*3600; % in seconds
%tend = 72*3600; %interval is shifted, see PD_strat
%tend = 12*3600;
tend = 84;%*3600; % in s

% Initial values

% y0 = [9.906e1; 6.624e8; 3*5.326e11; 2*1.697e16; 4.000e6; 2*1.093e9];
y0 = [9.906e+1; 6.624e+8; 5.326e+11; 1.697e+16; 4.000e+6; 1.093e+9];
n1 = [1;1;3;2;1;2];
y0 = y0 .* n1;

scaleode = 1e-8;


y0 = scaleode*y0;

% Time step size
dt0 = 3600;
dt0 = 0.01;
% dt0 = 4.5;

% Visualization amplification factors
%visfac = [1e-2 1e-8 1e-11 1e-11 1e-7 1e-9];
%visfac = [1e-2 1e-8 1e-11 1e-16 1e-16 1e-12];

% Plot range
%yrange =  [-1 11];

% Compound composition matrix
%E = [1 0 0 1 1 1; 0 1 1 0 0 0; 1 1 1 1 1 1];
E = [1 1 1 1 1 1];
E = [];

% Relaxation
entropy_flag = true;
if entropy_flag
    % n2 = [0;0;0;0;1;1];
    % if exist('n1','var')
    %     n2 = n2 ./ n1;
    % end
        % Entropy
    eta = @(Y) (Y(5,:) + 0.5*Y(6,:));
    % para = scaleode;

    entropy_cons = true;
    % Entropy prime
    eta_prime = @(y)[0;0;0;0;1;0.5];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 04/23/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, Y] = ode45wrapper(PD, y0, dt0, tint, varargin)

if nargin > 4,
  RTOL = varargin{1}(1); ATOL = varargin{1}(2);
else
  RTOL = 1e-2; ATOL = 1e-2;
end

if nargin > 5 && ~isa(varargin{2},'function_handle')
  nonneg = varargin{2};
else
  nonneg = [];
end

[~,~,f] = PD(0,y0);

fprintf('%s\n','ODE45')

options.InitialStep = dt0;
options.RelTol = RTOL;
options.AbsTol = ATOL;
options.NonNegative = nonneg;
options.Stats = 'on';
[t,Y] = ode45(f,tint,y0,options);
Y = Y';
t = t';

fprintf('\n')
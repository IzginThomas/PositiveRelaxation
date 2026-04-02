%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 04/23/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, Y] = ode23swrapper(PD, y0, dt0, tint, varargin)

if nargin > 4,
  RTOL = varargin{1}(1); ATOL = varargin{1}(2);
else
  RTOL = 1e-2; ATOL = 1e-2;
end

[~,~,f] = PD(0,y0);

fprintf('%s\n','ODE23s')

options.InitialStep = dt0;
options.RelTol = RTOL;
options.AbsTol = ATOL;
options.Stats = 'on';
[t,Y] = ode23s(f,tint,y0,options);
Y = Y';
t = t';

fprintf('\n')

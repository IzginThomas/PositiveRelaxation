%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)             %%%
%%% Date: 1/20/2025                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DENSE OUTPUT YET TO BE IMPLEMTED!
function [t, Y,str,varargout] = LW_linad(PD, y0, dt0, tint, varargin)
warning('off','all')
% Optional arguments
getpara();
order = 2; % expected order in time (and space) of the scheme

setdefaults(); % adap defaults are set in getpara.m

%get string for plot
str = sprintf('%s','LW2');


tstart = tint(1);
tend = tint(2);

t = tstart;
dt = dt0;
y = y0;
Y = y0;

cnt_rej = 0;
lengtht = 0;

while t(end) < tend

    if t(end) + dt > tend
        dt = tend - t(end);
    end
    yn = y;
    lambda = para*dt/dx;

    yn_m = circshift(yn, 1);   % u_{j-1}
    yn_p = circshift(yn, -1);  % u_{j+1}

    y = yn - 0.5 * lambda * (yn_p - yn_m) + ...
                0.5 * lambda^2 * (yn_p - 2*yn + yn_m);

    t(end+1) = t(end) + dt;
    Y(:,end+1) = y;


    % if length(t) == 1e+6  || dt < 1e-100
    %     lengtht = NaN;
    %     break
    % else
    %     lengtht = 0;
    % end

    fprintf('Current time in percent: %f\n',100*(t(end)-t(1))/(tend-tstart))
end
if nargout > 1
    varargout{1} = cnt_rej;
    varargout{2} = lengtht;
    varargout{3} = [];
    varargout{4} = order; % output for convergence test
end
%fprintf('%s(%.3f)\n',mfilename,a21)
%fprintf('%i successful steps\n',length(t)-1)
%fprintf('%i failed attempts\n\n',cnt_rej);

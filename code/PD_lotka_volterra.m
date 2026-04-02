%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)             %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, D, f, S, r, Rrem,Rminus] = PD_lotka_volterra(t,y,varargin)
getpara();


D(1,2) = y(1) * y(2);
P = D';



if nargout > 2
    f = @(t,y) [2*y(1) - y(1)*y(2); y(1)*y(2) - y(2)];
end

if nargout > 3 %% Need to be done
 S=zeros(length(y));
  
  
   r = @(y) zeros(length(y),1);

end

if nargout > 5 
Rrem = zeros(length(y),1); % remainder
Rminus = Rrem;
Rrem(1) = 2*y(1);
Rminus(2) = y(2);
end


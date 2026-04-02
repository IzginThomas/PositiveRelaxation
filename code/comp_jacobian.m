%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 09/02/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Df = comp_jacobian(f,t,y)

h = 1e-8;
N = length(y);

Df = zeros(N);
for i=1:N
  e = zeros(N,1);
  e(i) = 1;
  Df(:,i)=(f(t,y + h*e) - f(t,y - h*e))/(2*h);
  %Df(:,i)=(f(t,y + h*e) - f(t,y))/h;
end

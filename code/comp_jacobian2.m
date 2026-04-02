%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 11/03/2019                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Df = comp_jacobian2(F,x)

h = 1e-8;
N = length(x);

Df = zeros(N);
for i=1:N
  e = zeros(N,1);
  e(i) = 1;
  %Df(:,i)=(F(x + h*e) - F(x - h*e))/(2*h);
  Df(:,i)=(F(x + h*e) - F(x))/h;
end

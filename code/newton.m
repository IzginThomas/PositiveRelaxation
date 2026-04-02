%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 03/11/2019                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,cnt] = newton(F,x0,tol,maxit)
x = x0;

Fx = F(x);
res = norm(Fx);

cnt = 0;
while res > tol && cnt < maxit
  DF = comp_jacobian2(F,x);
  
  % update
  x = x + DF\-Fx;
  
  Fx = F(x);
  res = norm(Fx);
  
  cnt = cnt + 1;
end
if cnt == maxit && res > tol
  warning('ODESolverTestSuite:newton:Maximum number of nonlinear iterations reached!\n %5i %e %e\n',cnt,res)
end
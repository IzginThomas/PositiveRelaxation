%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 08/29/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = comp_err(t,Y,Yref,p,idx,relerr)
% N = size(Y,1);
%
% err = 0;
% for j = 1:N
%   err = err + sqrt(sum(diff(t).*abs(Yref(2:end,j)' - Y(j,2:end)).^2))/sqrt(sum(diff(t)'.*abs(Yref(2:end,j)).^2));
%   %err = err + sqrt(sum(diff(t).*abs(Yref(2:end,j)' - Y(j,2:end)).^2));
% end
% err = err/N;

Yref = Yref';

if p ~= inf
  aux = sum(abs(Y(idx,:) - Yref(idx,:)).^p,1);
  lpnorm_Ydiff = (sum(diff(t).*(0.5*(aux(1:end-1)+aux(2:end)))))^(1/p);
  err = lpnorm_Ydiff;
  
  if relerr
    aux = sum(abs(Yref(idx,:)).^p,1);
    lpnorm_Yex = (sum(diff(t).*(0.5*(aux(1:end-1)+aux(2:end)))))^(1/p);
    err = err/lpnorm_Yex;
  end
else
  linfnorm_Ydiff = max(max(abs(Y(idx,:) - Yref(idx,:))));
  err = linfnorm_Ydiff;
  
  if relerr
    linfnorm_Yex = max(max(abs(Yref(idx,:))));    
    err = err/linfnorm_Yex;
  end
end







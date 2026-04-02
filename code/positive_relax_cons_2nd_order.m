function [gamma, y_gamma, flag_relax] = positive_relax_cons_2nd_order(y0, yn, gamma, y_gamma, y_func, M_func, relax_tol, cnt_rej_relax,varargin)
getpara();
eta_ygamma = eta(y_gamma);
cnt = 0;
etaold = eta(y0); % better than eta(yn) due to error accumulation
while max(abs(eta_ygamma - etaold)) > relax_tol && gamma > 0 && ~isinf(gamma)
    [gamma, y_gamma, eta_ygamma] = positive_relax_cons_2nd_order_step(y_gamma,yn, gamma, eta, eta_ygamma,etaold,eta_prime,y_func,M_func);
     cnt = cnt + 1;
     if cnt == 1000
         break;
     end
end
if ((abs(eta_ygamma - etaold) <= relax_tol && ~isinf(gamma)) || cnt_rej_relax >= 1e3) && gamma > 0
   flag_relax = true;
else
    flag_relax = false;
end

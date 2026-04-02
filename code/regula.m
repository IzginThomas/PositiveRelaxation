function [gamma,flag_relax] = regula(r,lgamma,rgamma,rlgamma,rrgamma, relax_tol,max_it)

gamma1 = (lgamma * rrgamma - rgamma * rlgamma)/(rrgamma - rlgamma);

rgamma1 = r(gamma1);
i = 0;
while abs(rgamma1) > relax_tol  && i <= max_it && abs(rgamma-lgamma) > relax_tol
    if rgamma1 * rlgamma > 0
        lgamma = gamma1;
        rlgamma = rgamma1;
    elseif rgamma1 * rrgamma > 0
        rgamma = gamma1;
        rrgamma = rgamma1;
    else
        error('Requirments for convergence of Regular falsi not met')
    end
    gamma1 = (lgamma * rrgamma - rgamma * rlgamma)/(rrgamma - rlgamma);
    rgamma1 = r(gamma1);
    i = i +1;
end
gamma = gamma1;

if ((abs(rgamma1) <= relax_tol && ~isinf(gamma))) && gamma > 0
    flag_relax = true;
else
    flag_relax = false;
end
end
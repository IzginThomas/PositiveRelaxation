function [gamma, y_gamma,sigmagamma,flag_relax] = positive_relax_cons_secant(etaygammaold, gammaold, gamma, sigma_func, eta, eta_ygamma,etaold,y_func,relax_tol,max_it)
cnt = 0;
while abs(eta_ygamma - etaold) > relax_tol && gamma > 0 && ~isinf(gamma) && abs(gamma-gammaold)>relax_tol && cnt <= max_it
    gamma1 = gamma - (eta_ygamma - etaold)*(gamma - gammaold) / (eta_ygamma - etaygammaold);

    gammaold = gamma;
    etaygammaold = eta_ygamma;

    gamma = gamma1;
    if gamma > 0
        sigmagamma = sigma_func(gamma);
        y_gamma = y_func(gamma,sigmagamma);
        eta_ygamma = eta(y_gamma);
    else
        sigmagamma = 0;
        y_gamma = 0;
    end
    cnt = cnt + 1;
end
if abs(eta_ygamma - etaold) <= relax_tol && ~isinf(gamma) && gamma > 0
    flag_relax = true;
else
    flag_relax = false;
end

end
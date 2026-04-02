function  [gamma, y_gamma, sigmagamma, flag_relax] = relax_newton3(y_gamma,yn,sigmagamma, dmu_func,gamma, sigma_func, ...
                         eta, eta_ygamma,etaold,eta_prime,y_func,M_func, M4_func,mu_func,relax_tol, max_it)
    
gammaold = 0;
cnt = 0;
while abs(eta_ygamma - etaold) > relax_tol && gamma > 0 && ~isinf(gamma) && abs(gamma-gammaold)>relax_tol && cnt <= max_it
gammaold = gamma;
cnt = cnt + 1;
dsigmagamma = get_derivative(yn,sigmagamma,gamma,dmu_func(gamma),mu_func(gamma),M4_func(gamma));

dygamma = get_derivative(yn, y_gamma, gamma,dsigmagamma,sigmagamma,M_func(gamma,sigmagamma));

rPrime = (eta_prime(y_gamma))' * dygamma;

gamma = gamma - (eta_ygamma - etaold) / (rPrime + realmin);
if gamma > 0
    sigmagamma = sigma_func(gamma);
    y_gamma = y_func(gamma,sigmagamma);
    eta_ygamma = eta(y_gamma);
end

end
if abs(eta_ygamma - etaold) <= relax_tol && ~isinf(gamma) && gamma > 0
    flag_relax = true;
else
    flag_relax = false;
end

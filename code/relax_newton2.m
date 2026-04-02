function  [gamma, y_gamma, sigmagamma, flag_relax] = relax_newton2(y_gamma,yn,sigmagamma,gamma, sigma_func,dsigma_func, ...
                         eta, eta_ygamma,etaold,eta_prime,y_func,M_func,relax_tol, max_it,pos_ind)
ynew = y_gamma;
gammaold = 0;
cnt = 0;
id_exp = setdiff(1:length(yn),pos_ind);
while abs(eta_ygamma - etaold) > relax_tol && gamma > 0 && ~isinf(gamma) && abs(gamma-gammaold)>relax_tol && cnt <= max_it
gammaold = gamma;
cnt = cnt + 1;
dsigmagamma = dsigma_func(gamma);
dygamma = zeros(size(yn));
dygamma(pos_ind) = get_derivative(yn(pos_ind), y_gamma(pos_ind), gamma,dsigmagamma,sigmagamma,M_func(gamma,sigmagamma));
dygamma(id_exp) = ynew(id_exp) - yn(id_exp);

rPrime = (eta_prime(y_gamma))' * dygamma;

gamma = gamma - (eta_ygamma - etaold) / (rPrime + realmin);
if gamma > 0
    sigmagamma = sigma_func(gamma);
    y_gamma(pos_ind) = y_func(gamma,sigmagamma);
    y_gamma(id_exp) = yn(id_exp) + gamma *(ynew(id_exp) - yn(id_exp));
    eta_ygamma = eta(y_gamma);
end

end
if abs(eta_ygamma - etaold) <= relax_tol && ~isinf(gamma) && gamma > 0
    flag_relax = true;
else
    flag_relax = false;
end

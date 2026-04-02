function [gamma, y_gamma, sigmagamma, flag_relax] = positive_relax_cons_2(y0, yn, gamma, y_gamma, y_func, M_func,sigmagamma, sigma_func,dsigma_func, relax_tol, cnt_rej_relax, varargin)
getpara();
flag_relax = true;
eta_ygamma = eta(y_gamma);
max_it = 100;
etaold = eta(y0); % better than eta(yn) due to error accumulation
r = @(Gamma) eta(y_func(Gamma,sigma_func(Gamma))) - etaold;
switch relax_meth
    case 'newton'
        [gamma, y_gamma, sigmagamma, flag_relax] = relax_newton2(y_gamma,yn,sigmagamma,gamma, sigma_func,dsigma_func, ...
            eta, eta_ygamma,etaold,eta_prime,y_func,M_func,relax_tol, max_it,pos_ind);
    case 'secant'
        gammaold = 1.1 * abs(gamma);
        etaygammaold = eta(y_func(gammaold,sigma_func(gammaold)));
        [gamma, y_gamma, sigmagamma, flag_relax] = positive_relax_cons_secant(etaygammaold, gammaold, gamma, sigma_func, eta, eta_ygamma,etaold,y_func,relax_tol,max_it);
    otherwise
        if strcmp(relax_meth,'bisection') || strcmp(relax_meth,'regula')
            lgamma = 0.99 * gamma;
            rgamma = 1.01 * gamma;
            cnt_lr = 0;
            rlgamma = r(lgamma);
            rrgamma = r(rgamma);
            while rlgamma * rrgamma > 0 && abs(rrgamma) >= relax_tol && abs(rlgamma) >= relax_tol
                lgamma = 0.99 * lgamma;
                rgamma = 1.01 * rgamma;
                rlgamma = r(lgamma);
                rrgamma = r(rgamma);
                cnt_lr = cnt_lr + 1;
                if cnt_lr > max_it & rlgamma * rrgamma > 0
                    fprintf('No suitable bounds found after %d iterations.\n Left = %e; r(left) = %e \n Right = %e; r(right) = %e \n',...
                        cnt_lr, lgamma, rlgamma, rgamma, rrgamma);
                    flag_relax = false;
                    break;
                end
            end
            if flag_relax
                if strcmp(relax_meth,'bisection')
                    [gamma,flag_relax] = bisectionMethod(r, lgamma, rgamma, rlgamma, relax_tol,max_it);
                else
                    if abs(rrgamma) < relax_tol || abs(rlgamma) < relax_tol
                        if abs(1-lgamma) <= abs(1-rgamma)
                            gamma = lgamma;
                        else 
                            gamma = rgamma;
                        end
                        flag_relax = true;
                    else
                        [gamma,flag_relax] = regula(r,lgamma,rgamma,rlgamma,rrgamma,relax_tol,max_it);
                    end
                end
                sigmagamma = sigma_func(gamma);
                y_gamma = y_func(gamma,sigmagamma);
            end
        else
            error('Unkown relaxation method solver')
        end
end

end

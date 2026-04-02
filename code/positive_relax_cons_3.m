function [gamma, y_gamma, sigmagamma, flag_relax] = positive_relax_cons_3(y0, yn, sigmagamma, ...
    gamma, y_gamma, sigma_func, dmu_func, y_func, M_func, M4_func, mu_func,...
    relax_tol, varargin)
getpara();
flag_relax = true;
eta_ygamma = eta(y_gamma);
etaold = eta(y0); % better than eta(yn) due to error accumulation
max_it = 100;
r = @(Gamma) eta(y_func(Gamma,sigma_func(Gamma))) - etaold;
switch relax_meth
    case 'newton'
        [gamma, y_gamma, sigmagamma, flag_relax] = relax_newton3(y_gamma,yn,sigmagamma, dmu_func,gamma, sigma_func, ...
            eta, eta_ygamma,etaold,eta_prime,y_func,M_func, M4_func,mu_func,relax_tol,max_it);
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


%  gammaold = 0;
%  etaygammaold = eta(yn);
% while abs(eta_ygamma - etaold) > relax_tol && gamma > 0 && ~isinf(gamma)
%      %
%
%
%
%      cnt = cnt + 1;
%      if cnt == 100
%          break;
%      end
% end

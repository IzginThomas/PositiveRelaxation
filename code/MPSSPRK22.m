%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 04/02/2026                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, Y,str,varargout] = MPSSPRK22(PD, y0, dt0, tint, varargin)
getpara();
order = 2;
Gamma_vec = [];
a21 = parameter(1);
b10 = parameter(2);
str = sprintf('%s(%.2f,%.2f)','MPSSPRK2',a21,b10);

% a21 = alpha;
a10 = 1;
a20 = 1 - a21;

% b10 = beta;
b20 = 1 - 1/(2*b10) - a21*b10;
b21 = 1/(2*b10);

s = (b20 + b21 + a21*b10^2)/(b10*(b20+b21));

assert(a21 >= 0)
assert(a21 <= 1)
assert(b10 > 0)
assert(a21*b10 + 1/(2*b10) <= 1,sprintf('a21*b10 + 1/(2*b10) = %f > 1\n',a21*b10 + 1/(2*b10)))

% I = length(y0);
Idmat = speye(length(y0));

tstart = tint(1);
tend = tint(2);

t = tstart;
dt = dt0;
cnt_rej_relax = 0;

setdefaults(); % adap defaults are set in getpara.m
Y = y0;
y = y0 + realmin;
if relaxation
    relflag = true;
else
    relflag = false;
end
while t(end) < tend
    if t(end) + dt > tend
        dt = tend - t(end);
        relaxation = false;
    elseif relflag
        relaxation = true;
    end
    yn = y;
    [P,D,f,~,~,Rrem, Rminus] = PD(t(end),yn,varargin{:});
    M2 = Idmat + b10*dt*(diag(sum(D,2) + Rminus) - P)*diag(1./yn);
    y2 = M2\(a10*yn + b10*dt*Rrem);
    y2 = y2 + realmin;

    [P2, D2,~,~,~,R2rem,R2minus] = PD(t(end)+b10*dt,y2,varargin{:});

    sigma_func = @(gamma) y2.^(gamma*s).*yn.^(1-gamma*s) + realmin;
    dsigma_func = @(gamma) sigma_func(gamma).*s.*log(y2./yn);

    M_aux = (b20*dt*(diag(sum(D,2) + Rminus) - P) ...
        + b21*dt*(diag(sum(D2,2) + R2minus) - P2) );
    M = @(sigmagamma) M_aux*diag(1./sigmagamma);
    M_func = @(Gamma,sigmagamma) Gamma*M(sigmagamma) + Idmat;
    rem_func = @(Gamma) Gamma*dt*(b20*Rrem + b21*R2rem) + realmin;
    y_func = @(Gamma,sigmagamma) M_func(Gamma,sigmagamma) \ (yn + Gamma*a21*(y2-yn) + rem_func(Gamma)) + realmin;
    sigma = sigma_func(1) + realmin;
    y_aux = y_func(1,sigma);

    gamma = 1;
    if relaxation && entropy_flag
        relax_tol = 1e-13;
        entropy_tol = 1e-13;
        Gamma_vec(1) = gamma;
        if entropy_cons
            % Initialize
            y_gamma = y_aux;
            [gamma, y_gamma, ~, flag_relax] = positive_relax_cons_2(y0, yn, gamma, y_gamma, y_func, M_func,sigma, sigma_func,dsigma_func, relax_tol, cnt_rej_relax, varargin{:});
        else
            etaold = eta(yn);
            f1 = f(t(end), yn);
            f2 = f(t(end) + b10*dt, y2);
            b1 = a21 * b10 + b20;
            b2 = b21;
            aux_rhs = b1*(eta_prime(yn))'*f1 + b2*(eta_prime(y2))'*f2;

            etanew = etaold + dt* aux_rhs;

            ynew = y_aux;
            y_gamma = ynew; % gamma = 1
            eta_ygamma = eta(y_gamma);
            if eta_ygamma > etanew + entropy_tol % solve equation
                residual = abs(eta_ygamma - etaold - gamma*(etanew - etaold));
                % eta_prev = eta_ygamma;
                cnt = 0;
                while residual > relax_tol   && gamma > 0&& cnt <= 100
                    cnt = cnt + 1;
                    rPrime = (eta_prime(y_gamma))' * (ynew - yn) - (etanew - etaold);
                    if rPrime<-eps || (eta_prime(yn))' * (ynew - yn) - (etanew - etaold) > eps
                        % disp(dt);
                        gamma = -1;
                        break;
                    end
                    gamma = gamma - (eta_ygamma - etaold - gamma * (etanew - etaold)) / (rPrime + realmin);
                    y_gamma = yn + gamma * (ynew - yn);
                    eta_ygamma = eta(y_gamma);
                    residual = abs(eta_ygamma - etaold - gamma*(etanew - etaold));
                end
            end
            gamma = min(gamma,1);
            if (eta_ygamma <= etanew + entropy_tol || cnt_rej_relax >= 1e3) && gamma > 0
                flag_relax = true;
            else
                flag_relax = false;
            end
        end


        if flag_relax
            t(end + 1) = t(end) + gamma*dt;
            y = y_gamma;
            Y(:,end+1) = y;

                dt = 1.01 * dt;

            Gamma_vec(end + 1) = gamma;
            cnt_rej_relax = 0;
            flag = false; % for printing current progress

        else
            dt = dt * 0.9;
            % disp(dt);
            flag = true; % for not printing current progress

            cnt_rej_relax = cnt_rej_relax + 1;
        end
    else
        if ~exist('Gamma_vec','var')
            Gamma_vec = [];
        end
        t(end+1) = t(end) + dt;
        y = y_aux;
        Y(:,end+1) = y;
        flag = false; % for printing current progress
        %dt = CFL * dx;%/(max(para + abs(y(ind2))));

    end
    if ~flag
        fprintf('Current time in percent: %f\n',100*(t(end)-t(1))/(tend-tstart))
    end
end
if nargout > 1
    varargout{1} = 0;
    varargout{2} = length(t);
    varargout{3} = Gamma_vec;
    varargout{4} = order;
end
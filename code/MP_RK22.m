%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 1/20/2026                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DENSE OUTPUT YET TO BE IMPLEMTED!
function [t, Y,str,varargout] = MP_RK22(PD, y0, dt0, tint, varargin)
warning('off','all')
order = 2;
% Optional arguments
getpara();
a21 = parameter;

ind1 = pos_ind;
ind2 = setdiff(1:length(y0),ind1);

PWD_flag = 1; % for different dense output formulae

if PWD_flag
    q = @(gamma) 1/a21;
    dq = @(gamma) 0;
else
    q = @(gamma) gamma/a21;
    dq = @(gamma) 1/a21;
end


if exist('adap','var') && ~exist('Beta1','var')
    % C = {1.5193, -0.42576,-0.078535, -0.29465, 2}; % optimized controller
    C = {0.7, -0.4, 0, 0, 1};
    [Beta1, Beta2, Beta3, Alpha2, kappa2] = C{:};
end
setdefaults(); % adap defaults are set in getpara.m
skip_adap = false;
%get string for plot
if ~strcmp(errctrl,'off')
    if freq == 0
        str = sprintf('%s(%.3f)','MP_RK22adap',a21);
    else
        str = sprintf('%s(%.3f,%i)','MP_RK22adap',a21,freq);
    end
else
    if freq == 0
        str = sprintf('%s(%.3f)','MP_RK22',a21);
    else
        str = sprintf('%s(%.3f,%i)','MP_RK22',a21,freq);
    end
end


b2 = 1/(2*a21);
b1 = 1 - b2;

assert(a21 ~= 0)

tstart = tint(1);
tend = tint(2);

t = tstart;
dt = dt0;
dtold = dt;

Idmat = speye(length(ind1)); % ind1 = Modified Patankar-treated parts, ind2 = explicitly treated

if strcmp(errctrl,'PID') && t(end) == tstart
    epsminus = 1; %eps_{-1}
    eps0 = 1; %eps_0
end


y = y0 + realmin;
Y = y;

flag_reject_dt = true;
cnt_rej = 0;
cnt_rej_relax = 0;
lengtht = 0;
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
    y = y + realmin;
    yn = y;

    yn1 = yn(ind1);

    yn2 = yn(ind2);

    [P,D,f,~,~,Rrem, Rminus] = PD(t(end),yn,varargin{:});
    y2 = yn; % initialize
    M2 = Idmat + a21*dt*(diag(sum(D,2) + Rminus) - P)*diag(1./yn1);
    y2(ind1) = M2\(yn1 + a21*dt*Rrem) + realmin;

    fn = f(t(end),yn);
    fn2 = fn(ind2);
    y2(ind2) = yn2 + a21*dt*fn2;

    [P2, D2,~,~,~,R2rem,R2minus] = PD(t(end)+a21*dt,y2,varargin{:});

    sigma_func = @(gamma) yn1.*(y2(ind1)./yn1).^(q(gamma)) + realmin;
    dsigma_func = @(gamma) sigma_func(gamma).*dq(gamma).*log(y2(ind1)./yn1);

    sigma = sigma_func(1);
    y_aux = zeros(size(y0));

    M_aux = (b1*dt*(diag(sum(D,2) + Rminus) - P) ...
        + b2*dt*(diag(sum(D2,2) + R2minus) - P2) );
    M = @(sigmagamma) Idmat + M_aux*diag(1./sigmagamma);
    M_func = @(Gamma,sigmagamma) Gamma*(M(sigmagamma) - Idmat) + Idmat;
    rem_func = @(Gamma) Gamma*dt*(b1*Rrem + b2*R2rem) + realmin;
    y_func = @(Gamma,sigmagamma) M_func(Gamma,sigmagamma) \ (yn1 + rem_func(Gamma)) + realmin;
    y_aux(ind1) = y_func(1,sigma);

    f2 = f(t(end)+a21*dt,y2);
    f22 = f2(ind2);
    y_aux(ind2) = yn2 + dt*(b1*fn2 + b2*f22);

    gamma = 1;

    if relaxation && entropy_flag
        relax_tol = 1e-13;
        entropy_tol = 1e-13;
        Gamma_vec(1) = gamma;
        if entropy_cons
            % Initialize
            y_gamma = y_aux;

            [gamma, y_gamma, sigma, flag_relax] = positive_relax_cons_2(y0, yn, gamma, y_gamma, y_func, M_func,sigma, sigma_func,dsigma_func, relax_tol, cnt_rej_relax, varargin{:});
        else
            etaold = eta(yn);
            f1 = f(t(end), yn);
            f2 = f(t(end) + a21*dt, y2);
            aux_rhs = b1*(eta_prime(yn))'*f1 + b2*(eta_prime(y2))'*f2;

            etanew = etaold + dt* aux_rhs;

            ynew = y_aux;
            y_gamma = ynew; % gamma = 1
            sigmagamma = sigma; % gamma = 1
            eta_ygamma = eta(y_gamma);
            if eta_ygamma > etanew + entropy_tol % solve equation
                residual = abs(eta_ygamma - etaold - gamma*(etanew - etaold));
                cnt = 0;
                while residual > relax_tol   && gamma > 0 && cnt <= 100
                    cnt = cnt + 1;
                    rPrime = (eta_prime(y_gamma))' * (ynew - yn) - (etanew - etaold);

                    if rPrime <=-eps ||  (eta_prime(yn))' * (ynew - yn) - (etanew - etaold) > eps
                        % disp(norm(f1));
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
            if strcmp(errctrl,'off')
                % dt = gamma * dt;
                dt = 1.01 * dt;
                %dt = dt0;

                %dt = CFL * dx;%/(max(para + abs(y(ind2))));
            else
                % if ~entropy_cons
                sigma = sigma_func(gamma);
                % end
                skip_adap = false;
                flag_reject_dt = false;
            end
            Gamma_vec(end + 1) = gamma;
            cnt_rej_relax = 0;
            flag = false; % for printing current progress
        else
            dt = dt * 0.9;
            flag = true; % for not printing current progress
            if ~strcmp('errctrl','off')
                skip_adap = true;
            end
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

    if ~skip_adap
        if strcmp(errctrl,'PID')
            dt_temp = dt;
            [dt,flag,eps1] = compute_dt(y,yn,[sigma;yn2+a21*dt*fn2],dt,RTOL,ATOL,errctrl,order,epsminus,eps0,Beta1,Beta2,Beta3, Alpha2, kappa2, dtold);
            if ~flag
                dtold = dt_temp;
            end
        elseif strcmp(errctrl,'off')
            flag = false;
        else
            [dt,flag] = compute_dt(y,yn,sigma,dt,RTOL,ATOL,errctrl,order);
        end
        if flag & flag_reject_dt
            Y(:,end) = [];
            t(end) = [];
            y = Y(:,end);
            cnt_rej = cnt_rej + 1;
            if exist('Gamma_vec','var') && ~isempty(Gamma_vec)
                Gamma_vec(end) = [];
            end
            % disp(cnt_rej)
        else
            if strcmp(errctrl,'PID')
                epsminus = eps0;
                eps0 = eps1;
            end
        end
    else
        epsminus = 1;
        eps0 = 1;
        skip_adap = false;
    end
    if ~flag
        fprintf('Current time in percent: %f\n',100*(t(end)-t(1))/(tend-tstart))
    end
end
if nargout > 1
    varargout{1} = cnt_rej;
    varargout{2} = lengtht;
    varargout{3} = Gamma_vec;
    varargout{4} = order;
end
%fprintf('%s(%.3f)\n',mfilename,a21)
%fprintf('%i successful steps\n',length(t)-1)
%fprintf('%i failed attempts\n\n',cnt_rej);

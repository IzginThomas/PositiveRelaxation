%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 1/20/2026                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [t, Y,str,varargout] = MPRK43I(PD, y0, dt0, tint, varargin)
% warning('off','all')
getpara();
alpha = parameter(1);
beta = parameter(2);
small_const = realmin;%1e-8;
if exist('adap','var') && ~exist('Beta1','var')
    C = {1.7706, -0.27744, -0.37701, -0.95947, 3};
    [Beta1, Beta2, Beta3, Alpha2, kappa2] = C{:};
end
setdefaults(); % adap defaults are set in getpara.m
skip_adap = false;
%get string for plot
if ~strcmp(errctrl,'off')
    if freq == 0
        str = sprintf('%s(%.3f,%.3f)','MPRK43Iadap',alpha,beta);
    else
        str = sprintf('%s(%.3f,%.3f,%i)','MPRK43Iadap',alpha,beta,freq);
    end
else
    if freq == 0
        str = sprintf('%s(%.3f,%.3f)','MPRK43I',alpha,beta);
    else
        str = sprintf('%s(%.3f,%.3f,%i)','MPRK43I',alpha,beta,freq);
    end
end

TOL = 1e-14;
assert(alpha >= 0.5 && alpha ~= 2/3);

alpha0 = 1/6*(3 + (3-2*sqrt(2))^(1/3) + (3+2*sqrt(2))^(1/3));
if (1/3<= alpha && alpha < 2/3)
    assert((beta - 2/3) >= 0 && (beta-3*alpha*(1-alpha))<=TOL)
elseif (2/3 < alpha && alpha < alpha0)
    assert(3*alpha*(1-alpha) <= beta && beta <= 2/3)
else
    assert((3*alpha-2)/(6*alpha-3)<=beta && beta<= 2/3)
end

a21 = alpha;
a31 = (3*alpha*beta*(1-alpha)-beta^2)/(alpha*(2-3*alpha));
a32 = (beta*(beta-alpha))/(alpha*(2-3*alpha));
b1 = 1 + (2-3*(alpha+beta))/(6*alpha*beta);
b2 = (3*beta-2)/(6*alpha*(beta-alpha));
b3 = (2-3*alpha)/(6*beta*(beta-alpha));


beta2 = 1/(2*a21);
beta1 = 1-beta2;

% bvec = {[beta1,beta2], [b1,b2,b3]};

p = 1/(3*a21*(a31+a32)*b3);
q = @(gamma) gamma/a21;

assert(abs(b2*a21+b3*(a31+a32)-1/2)<1e-14)
assert(abs(b2*a21^2+b3*(a31+a32)^2-1/3)<1e-14)
assert(abs(a21*a32*b3-1/6)<1e-14)
assert(abs(beta1+beta2-1)<1e-14)
assert(abs(a21*beta2-1/2)<1e-14)

N = length(y0);
Idmat = speye(length(y0));

tstart = tint(1);
tend = tint(2);

t = tstart;
dt = dt0;
dtold = dt;
%W = [];
% gammaint = linspace(0, 1, freq + 2);
% gammaint = gammaint(2:end);
Y = y0;
y = y0;
y = y + small_const;
cnt = 0;
cnt_rej = 0;
% cnt_rej_relax = 0;
lengtht = 0;
tn = t;
if strcmp(errctrl,'PID') && t(end) == tstart
    epsminus = 1; %eps_{-1}
    eps0 = 1; %eps_0
end
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
    yn = y + small_const;
    [P{1},D{1},f,~,~,Rrem{1}, Rminus{1}] = PD(t(end),yn,varargin{:});

    M2 = Idmat + a21*dt*(diag( sum(D{1},2) + Rminus{1} ) - P{1})*diag(1./yn);
    y2 = M2\(yn + a21*dt*Rrem{1})  + small_const;

    [P{2}, D{2},~,~,~,Rrem{2},Rminus{2}] = PD(t(end) + alpha * dt, y2, varargin{:});

    rho = yn.^(1-p).*(y2).^p + small_const;
    M3 = Idmat + (a31*dt*(diag( sum(D{1},2) + Rminus{1} ) - P{1}) ...
        + a32*dt*(diag( sum(D{2},2) + Rminus{2} ) - P{2}))*diag(1./rho);
    y3 = M3\(yn + dt*(a31*Rrem{1} + a32*Rrem{2}) + small_const) + small_const;

    mu_func =@(gamma) (yn).^(1-q(gamma)) .* (y2).^(q(gamma)) + small_const;
    M4_func =@(gamma) Idmat + gamma*(beta1*dt*(diag( sum(D{1},2) + Rminus{1} ) - P{1}) ...
        + beta2*dt*(diag( sum(D{2},2) + Rminus{2} ) - P{2}))*diag(1./mu_func(gamma));

    sigma_func = @(gamma) M4_func(gamma)\(yn + gamma*dt*(beta1*Rrem{1} + beta2*Rrem{2})+ small_const)  + small_const;


    [P{3}, D{3},~,~,~,Rrem{3},Rminus{3}] = PD(t(end) + beta*dt,y3, varargin{:});
    M_func = @(gamma,sigmagamma) Idmat + gamma*(b1*dt*(diag( sum(D{1},2) + Rminus{1} ) - P{1}) ...
        + b2*dt*(diag( sum(D{2},2) + Rminus{2} ) - P{2}) ...
        + b3*dt*(diag( sum(D{3},2) + Rminus{3} ) - P{3}))*diag(1./sigmagamma) ;

    y_func = @(gamma,sigmagamma)  M_func(gamma,sigmagamma)\(yn + gamma*dt*(b1*Rrem{1} + b2*Rrem{2} + b3*Rrem{3}) + small_const)  + small_const;
    sigma = sigma_func(1);
    y_aux = y_func(1,sigma);
    gamma = 1;
    if relaxation && entropy_flag
        relax_tol = 1e-13;
        Gamma_vec(1) = gamma;
        % etaold = eta(y0);
        if entropy_cons
            % InitializeR2minus
            y_gamma = y_aux;
            dmu_func = @(gamma) (mu_func(gamma) ./ a21) .* log(y2 ./ yn);
            [gamma, y_gamma, sigma, flag_relax] = positive_relax_cons_3(y0, yn, sigma, ...
                gamma, y_gamma, sigma_func, dmu_func, y_func, M_func, M4_func, mu_func,...
                relax_tol, varargin{:});
        else
            etaold = eta(yn);
            f1 = f(t(end), yn);
            f2 = f(t(end) + alpha*dt, y2);
            f3 = f(t(end) + beta*dt, y3);
            etanew = etaold + dt* (b1*(eta_prime(yn))'*f1 + b2*(eta_prime(y2))'*f2 + b3*(eta_prime(y3))'*f3);
            ynew = y_aux;
            y_gamma = ynew; % gamma = 1
            eta_ygamma = eta(y_gamma);
            if eta_ygamma > etanew % solve equation
                residual = abs(eta_ygamma - etaold - gamma*(etanew - etaold))/max(abs(eta_ygamma),abs(etaold));
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
            if (eta_ygamma <= etanew) && gamma >0 %|| cnt_rej_relax == 1e3) %&& dt0*gamma > 1e-3
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
                dt = 1.01 * dt;
                % dt = dt0;
            else
                skip_adap = false;
            end
            Gamma_vec(end + 1) = gamma;
        else
            dt = dt * 0.9;
            if ~strcmp('errctrl','off')
                skip_adap = true;
            end
        end
    else
        if ~exist('Gamma_vec','var')
            Gamma_vec = [];
        end
        t(end+1) = t(end) + dt;
        y = y_aux;
        Y(:,end+1) = y;
    end

    pp = 3; % order of main method
    if ~skip_adap
        if strcmp(errctrl,'PID')
            dt_temp = dt;
            [dt,flag,eps1] = compute_dt(y,yn,sigma,dt,RTOL,ATOL,errctrl,pp,epsminus,eps0, Beta1,Beta2,Beta3, Alpha2, kappa2, dtold);
            if ~flag
                dtold = dt_temp;
            end
        elseif strcmp(errctrl,'off') % No adaptive time step
            flag = false;
        else
            [dt,flag] = compute_dt(y,yn,sigma,dt,RTOL,ATOL,errctrl,p);
        end
        if flag
            Y(:,end) = [];
            t(end) = [];
            y = Y(:,end);
            cnt_rej = cnt_rej + 1;
            if exist('Gamma_vec','var') && ~isempty(Gamma_vec)
                Gamma_vec(end) = [];
            end
        else
            if strcmp(errctrl,'PID')
                epsminus = eps0;
                eps0 = eps1;
            end
        end
    else
        epsminus = 1;
        eps0 = 1;
    end

    if cnt_rej == 1e+4 || cnt_rej > 1e2 * length(t)
        cnt_rej = NaN;
        break
    end
    fprintf('Current time in percent: %f\n',100*(t(end)-t(1))/(tend-tstart))
end

if nargout > 1
    % varargout{1} = str;
    varargout{1} = cnt_rej;
    varargout{2} = lengtht;
    varargout{3} = Gamma_vec;
    % varargout{3} = deltat;
end
%fprintf('%s(%.3f,%.3f)\n',mfilename,alpha,beta)
%fprintf('%i successful steps\n',length(t)-1)
%fprintf('%i failed attempts\n\n',cnt_rej);

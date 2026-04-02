function [x_res,xlines] = globalNM2Deufel(f,x0,lam_start)

if nargin < 3 | isempty(lam_start), lam_start = 1; end

TOL     = 10^(-8);
lammin  = 10^(-15);

k = 0;
xk = x0;
xlines = [xk];

while k<10
   %Abbruchsbedingung
   if norm(f(xk))<TOL
        x_res= xk;
        break;
   end
    
    %Pr�dikotschritt
    df = jacob(f,xk);
    if cond(df)>=10^10
        if 0.9*lam_start <10^-7
            x_res = false
            break;
        end
        [x_res,xlines] = globalNM2Deufel(f,x0,0.9*lam_start);
        break
    end
    dk = df\-f(xk);
    
    if k == 0
        lam_k       = lam_start;        
    else
        sig_k_u     = norm(f(xk))/norm(f(xkm))*sig_k_u;
        lam_k       = min(1,1/sig_k_u);
    end
    
    %Regularitytest
    if lam_k < lammin
        %disp('Fehler der Konvergenz')
        x_res = false;
        break;
    end
    %Berchnung vorl�ufiges x_{k+1}
    %2.
    xkp = xk + lam_k * dk;
    fkp = f(xkp);
    %3.
    theta_k = norm(fkp)/norm(f(xk));
    sig_k_u     = 2*norm(f(xkp)-(1-lam_k)*f(xk))/(lam_k^2*norm(f(xk)));
    
    while theta_k >= 1-lam_k*1/4
        lam_k = min(1/sig_k_u,1/2*lam_k);
         if lam_k < lammin
            %disp('Fehler der Konvergenz')
            x_res = false;
            break;
         end
        xkp = xk + lam_k * dk;
        fkp = f(xkp);
        theta_k = norm(fkp)/norm(f(xk));
        sig_k_u     = 2*norm(f(xkp)-(1-lam_k)*f(xk))/(lam_k^2*norm(f(xk)));     
    end
    lam_k_s = min(1,1/sig_k_u);
    while lam_k_s > 4*lam_k
        lam_k = lam_k_s;
        xkp = xk + lam_k * dk;
        fkp = f(xkp);
        theta_k = norm(fkp)/norm(f(xk));
        sig_k_u     = 2*norm(f(xkp)-(1-lam_k)*f(xk))/(lam_k^2*norm(f(xk)));
    
        while theta_k >= 1-lam_k*1/4
            lam_k = min(1/sig_k_u,1/2*lam_k);
            if lam_k < lammin
                disp('Fehler der Konvergenz')
            break;
            end
        xkp = xk + lam_k * dk;
        fkp = f(xkp);
        theta_k = norm(fkp)/norm(f(xk));
        sig_k_u     = 2*norm(f(xkp)-(1-lam_k)*f(xk))/(lam_k^2*norm(f(xk)));     
        end
        lam_k_s = min(1,1/sig_k_u);   
    end
    
    
    xkm     = xk;
    xk      = xkp;
    k       = k+1;
    xlines   = [xlines,xk];
    
end   

% eingebundene funktionen
    function [J]=jacob(func,x)
    % computes the Jacobian of a function
    n=length(x);
    fx=feval(func,x);
    eps=1.e-8; % could be made better
    xperturb=x;
    for i=1:n
        xperturb(i)=xperturb(i)+eps;
        J(:,i)=(feval(func,xperturb)-fx)/eps;
        xperturb(i)=x(i);
    end
    end
  
end
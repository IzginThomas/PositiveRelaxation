if ~exist('dtfunc','var')
    dtfunc = @(dt) dt;
end
if ~exist('scaleode','var')
    scaleode = 1;
end
if ~exist('para','var')
    para = NaN;
end
if ~exist('d','var')
    d = 1;
end
if ~exist('relaxation','var')
   relaxation = false;
end
if ~exist('freq','var')
    freq = 0;
end
if ~exist('pde_plot','var')
    pde_plot = false;
end
if ~exist('eta','var') || ~exist('eta_prime','var')
    entropy_flag = false;
else
    entropy_flag = true;
    if ~exist('entropy_cons', 'var')
        entropy_cons = false;
    end
end
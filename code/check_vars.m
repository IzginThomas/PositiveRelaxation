if exist('dx','var')
    input = [input(:)', 'dx', dx];
end
if exist('Nx','var')
    input = [input(:)', 'Nx', Nx];
end
if exist('CFL','var')
    input = [input(:)', 'CFL', CFL];
end
if ~exist('pos_ind','var')
    pos_ind = 1:length(y0);
end
if ~exist('mass_cons_relax','var')
    mass_cons_relax = true;
end
input = [input(:)', 'mass_cons_relax', mass_cons_relax];

if ~exist('inds','var')
    inds = 1:length(y0);
end
input = [input(:)', 'pos_ind', pos_ind];
if exist('bc','var')
    input = [input(:)', 'bc', bc];
end
if exist('phys_flux','var')
    input = [input(:)', 'phys_flux', func2str(phys_flux)];
end
if exist('dphys_flux','var')
    input = [input(:)', 'dphys_flux', func2str(dphys_flux)];
end
if exist('Phi','var')
    input = [input(:)', 'Phi', func2str(Phi)];
end
if exist('PhiR','var')
    input = [input(:)', 'PhiR', func2str(PhiR)];
end
if exist('num_flux','var')
    input = [input(:)', 'num_flux', func2str(num_flux)];
end

if exist('scaleode','var')
    input = [input(:)', 'scaleode', scaleode];
else
    scaleode = 1;
end
if ~exist('dtfunc','var')
    dtfunc = @(dt) dt;
end
input = [input(:)', 'dtfunc', func2str(dtfunc)];
if exist('para','var')
    input = [input(:)', 'para', para];
end
if exist('jac_flux','var')
    input = [input(:)', 'jac_flux', func2str(jac_flux)];
end
if ~exist('entropy_flag','var')
    entropy_flag = false;
end
if entropy_flag
    input = [input(:)', 'entropy', func2str(eta),...
        'entropy_cons', entropy_cons,...
        'entropy_prime', func2str(eta_prime)];
end
if ~exist('pde_plot','var')
    pde_plot = false;
end
if ~exist('xscale','var')
    xscale = 'linear';
end
if ~exist('yscale','var')
    yscale = 'linear';
end
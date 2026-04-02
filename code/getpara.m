nInputs = numel(varargin);
if nInputs > 0
    for i = 1 : 2 : length(varargin{:})
        name = varargin{:}{i};
        value = varargin{:}{i+1};
        switch name
            case 'adap', adap = value;
            case 'dense', freq = value;
            case 'dtfunc', dtfunc = str2func(value);
            case 'dx', dx = value;
            case 'CFL', CFL = value;
            case 'pos_ind', pos_ind = value;
            case 'Nx', Nx = value;
            case 'bc', bc = value;
            case 'num_flux', num_flux = str2func(value);
            case 'jac_flux', jac_flux = str2func(value);
            case 'Phi', Phi = str2func(value);
            case 'PhiR', PhiR = str2func(value);
            case 'dphys_flux', dphys_flux = str2func(value);
            case 'phys_flux', phys_flux = str2func(value);
            case 'parameter', parameter = value;
            case 'para', para = value;
            case 'scaleode', scaleode = value;
            case 'relaxation', relaxation = value{1}; relax_meth = value{2};
            case 'entropy', eta = str2func(value);
            case 'entropy_prime', eta_prime = str2func(value);
            case 'entropy_cons', entropy_cons = value;
            case 'mass_cons_relax', mass_cons_relax = value;
            otherwise, error(['Invalid property:',name]);
        end
    end
    if exist('para','var') &&  exist('phys_flux','var') && nargin(phys_flux) > 1
        phys_flux = @(y)phys_flux(y,para);
    end
    if exist('para','var') &&  exist('dphys_flux','var') && nargin(dphys_flux) > 1
        dphys_flux = @(y) dphys_flux(y,para);
    end
    if exist('para','var') &&  exist('jac_flux','var') && nargin(jac_flux) > 1
        jac_flux = @(y)jac_flux(y,para);
    end
    if exist('para','var') &&   exist('eta','var') && nargin(eta) > 1
        eta = @(Y)eta(Y,para);
        eta_prime = @(y) eta_prime(y,para);
    end
    if exist('adap','var')
        if strcmp(class(adap),'cell')
            ladap = length(adap);
            if ladap == 1
                [RTOL, ATOL] = adap{1}{:};
            else
                if ~ max(strcmp(adap{1},'opt'))
                    [Beta1,Beta2,Beta3,Alpha2,kappa2]=adap{1}{:};
                end
                [RTOL, ATOL] = adap{2}{:};
                if ladap == 2
                    errctrl = 'PID';
                else % ladap >= 3
                    errctrl = adap{3};
                end
            end
        elseif isa(class(adap), 'char')
            RTOL = 1e-2;
            ATOL = 1e-2;
            errctrl = adap;
        end
        if ~exist('errctrl','var')
            errctrl = 'PID';
        end
    else
        errctrl = 'off';
    end
end
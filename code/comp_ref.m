
if ~exist('exsol', 'var')
    if ~isOctave
        options.AbsTol = 1e-13;
        options.RelTol = 1e-13;
        [~,~,f] = PD(tstart,y0,input);
        if strcmp(test,'mseir')
            %[tex,yex] = ode23s(f,tint,y0,options);
            [tex,yex] = ode15s(f,tint,y0,options);
        elseif strcmp(test,'bertolazzi')
            opt.RelTol = 1e-5;
            opt.AbsTol = opt.RelTol;
            opt.NonNegtaive = [0 1 1];
            [tex,yex] = ode23s(f,tint,y0,opt);
        else
            options.NonNegative = 1:length(y0);
            options.NonNegative = [];
            [tex,yex] = ode15s(f,tint,y0,options);
            %[tex,yex] = ode23s(f,tint,y0,options);
        end
    else
        warning('Octave is currently not supported.')
        return
    end
else
    if ~pde_plot
        tex = linspace(tstart,tend,(tend-tstart)*1000);
        yex = exsol(tex');
    else %PDE with exact solution is plotted at tend
        tex = x;
        yex = exsol(tend,x')';
    end
end
yex = yex/scaleode;

% linear invariants
if exist('E','var')
else
    E = ones(1,size(yex,2));
end

if isempty(E)
    sumyex = [];
else
    sumyex = yex*E';
end

yex = [yex sumyex];
if exist('visfac','var')
    yex = bsxfun(@times, yex, visfac);
end
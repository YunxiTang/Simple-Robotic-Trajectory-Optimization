function [x, res, H_f, free] = boxQP(H, g, lb, ub, x0, options)
%BOXQP A QP solver to solve box constraint
%      Minimize 0.5*x'*H*x + g'*x
%               s.t.
%               lb <= x <= ub
%   Inputs:
%       H : PD matrix (n x n)
%       g : bias vector (n)
%       lb: lower bound (n)
%       ub: upper bound (n)
%       x0: initial guess (n)
%       options: option for solver
%   Outputs:
%       x: solution (n)
%       result: result type (roughly, higher is better, see below)
%       H_f: subspace cholesky factor   (n_free * n_free)
%       free: set of free dimensions    (n)
    n = size(H,1);
    clamped = false(n,1);
    oldvalue = 0.0;
    res = 0;
    g_norm = 0.0;
    n_factor = 0.0;
    H_f = zeros(n);
    parser = struct();
    
    %% some functions
    clamp = @(x) max(lb, min(ub, x));
    obj_val = @(x) x'*g + 0.5*x'*H*x;
    Grad = @(x) g + H*x;
    
    %% set starting point
    if nargin > 4 && numel(x0)==n
        x = clamp(x0(:));
    else
        x = (lb + ub) / 2;
    end
    x(isnan(x)) = 0.0;
    
    %% set options
    if nargin > 5
        options = num2cell(options(:));
        [parser.maxIter, parser.minGrad, parser.minRelImprove, parser.stepDec, parser.minStep, parser.Armijo, parser.print] = deal(options{:});
    else % defaults
        parser.maxIter        = 100;       % maximum number of iterations
        parser.minGrad        = 1e-8;      % minimum norm of non-fixed gradient
        parser.minRelImprove  = 1e-8;      % minimum relative improvement
        parser.stepDec        = 0.8;       % factor for decreasing stepsize
        parser.minStep        = 1e-22;     % minimal stepsize for linesearch
        parser.Armijo         = 0.1;   	% Armijo parameter (fraction of linear improvement required)
        parser.print          = 0;			% verbosity
    end
    
    %% initial objective function
    value = obj_val(x);
    
    if parser.print > 0
        fprintf('[INFO]: Starting box-QP, dimension %-3d, initial value: %-12.3f\n',n, value);
    end
    
    %% Go-to main Loop
    for iter = 1:parser.maxIter
        if res ~= 0
            break
        end
        
        % check relative improvement
        if iter > 1 && (oldvalue - value) < parser.minRelImprove * abs(oldvalue)
            res = 4;
            break
        end
        oldvalue = value;
        
        % get gradient
        grad = Grad(x);
        
        % find clamped dimentions (constrained dimensions)
        old_clamped = clamped;
        clamped = false(n,1);              % logical variable
        clamped((x-lb)<1e-15 & (grad>0)) = true;
        clamped((ub-x)<1e-15 & (grad<0)) = true;
        free = ~clamped;
        
        % check if all clamped
        if all(clamped)
            res = 6;
            break
        end
        
        % factorize if clamped has changed
        if iter == 1
            factorize = true;
        else
            factorize = any(old_clamped ~= clamped);
        end
        
        if factorize
            [H_f, indef] = chol(H(free, free));
            if indef
                res = -1;
                break
            end
            n_factor = n_factor + 1;
        end
        
        % check gradient norm
        g_norm = norm(grad(free));
        if g_norm < parser.minGrad
            res = 5;
            break
        end
        
        % get search direction
        grad_clamped = Grad(x.*clamped);
        search = zeros(n,1);
        search(free)   = -H_f\(H_f'\grad_clamped(free)) - x(free);
        
        % check for descent direction
        sdotg = sum(search.*grad);
        if sdotg >= 0 % (should not happen)
            break
        end
        
        % armijo linesearch
        step  = 1;
        nstep = 0;
        
        xc = clamp(x+step*search);
        vc = obj_val(xc);
        
        while (vc - oldvalue) / (step*sdotg) < parser.Armijo
            step  = step*parser.stepDec;
            nstep = nstep+1;
            xc    = clamp(x+step*search);
            vc    = obj_val(xc);
            if step < parser.minStep
                res = 2;
                break
            end
        end
        
        if parser.print > 1
            fprintf('iter %-3d  value % -9.5g |g| %-9.3g  reduction %-9.3g  linesearch %g^%-2d  n_clamped %d\n', ...
                iter, vc, g_norm, oldvalue-vc, parser.stepDec, nstep, sum(clamped));
        end
        
        % accept candidate
        x     = xc;
        value = vc;
    end
    if iter >= parser.maxIter
        res = 1;
    end
end


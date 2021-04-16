function demoQP()
    options = [100 1e-8 1e-8 0.6 1e-22 0.1 2]; % defaults with detailed printing
    n 		= 15;
    rng(10);
    g 		= randn(n,1);
    H 		= randn(n,n);
    H 		= H*H';
    lower 	= -ones(n,1);
    upper 	=  ones(n,1);
    tic
    [x, res, H_f, free] = boxQP(H, g, lower, upper, randn(n,1), options);
    disp(x)
    toc
end
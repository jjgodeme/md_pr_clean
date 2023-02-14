function  z=perform_MDBT(z,gradf,f,gradpsi,psi,options)
% We implement here the backtracking version of MD.

options.null = 0;
niter = getoptions(options, 'niter', 100);
Lbnd  = getoptions(options,'Lbnd',3);
kappa = getoptions(options, 'kappa', 0);
eta   = getoptions(options, 'eta', 1.4);

for t=1:niter
  zprev=z;
  fprev=f(zprev);
  gradfprev=gradf(zprev);
  psiprev=psi(zprev);
  gradpsiprev=gradpsi(zprev);
  L = Lbnd;
  for i=1:1E2
    gamma=(1-kappa)/L;
    pgamma=gradpsiprev-gamma.*gradfprev;		% Computation of the temporary variable pgamma
    racine=hgroots([norm(pgamma,'fro')^2,0, 1,-1]);     % Find roots of the 3rd polynom   % Solvy 3rd degree polynom.
    z=racine*pgamma;
    Dpsi=psi(z)-psiprev-gradpsiprev'*(z-zprev);
    Df=f(z)-fprev-gradfprev'*(z-zprev);
    if Df>L*Dpsi
        break;
    else
        L=L/eta;
    end
  end
  gamma=(1-kappa)/L;
  pgamma=gradpsiprev-gamma.*gradfprev;			% Computation of the temporary variable pgamma
  racine=hgroots([norm(pgamma,'fro')^2,0, 1,-1]);	% Find roots of the 3rd polynom   % Solvy 3rd degree polynom;
  z=racine*pgamma;
end

function  z=perform_MD(z,gradf,gradpsi,options)
% We implement here standard Mirror Descent appproach

options.null = 0;
niter = getoptions(options, 'niter', 100);
kappa = getoptions(options, 'kappa', 0);
L     = getoptions(options,'L',3);
gamma = (1-kappa)/L;
for  t=1:niter
    pgamma=gradpsi(z)-gamma.*gradf(z);		                    % Computation of the temporary variable pgamma
    racine=hgroots([norm(pgamma,'fro')^2,0, 1,-1]);     % Find roots of the 3rd polynom   %Solvy 3rd degree polynom
    z     =racine*pgamma;
end
end

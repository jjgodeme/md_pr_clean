function z=perform_SG(z,gradf,f,niter)
% Perform the subgradients methods

  for t=1:niter
    fz=f(z);
    grad=gradf(z);
    d=fz*grad/norm(grad,'fro')^2;
    d(grad == 0) = 0;
    z=z-d;
   end 
end

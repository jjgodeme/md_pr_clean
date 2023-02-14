function z=perform_WF(z,gradf,options)
% Computes the gradient descent for PR

  tau0 = 330;                         	       % Time constant for step size
  mu   = @(t) min(1-exp(-t/tau0), 0.2);          % Schedule for step size
  niter   = getoptions(options, 'niter', 100);
  normest = getoptions(options, 'normest','error');

  for t=1:niter
    grad =gradf(z);
    z= z - mu(t)/normest^2 * grad;
  end
end

function [zinit,normest]= Spectral_init(z0,Y,A,At,npower_iter)
% The following function computes the spectral initialization of the
% Phase Retrieval problem using a random guess z0
  z0= z0/norm(z0,'fro');
  for tt = 1:npower_iter
    z0 = At(Y.*A(z0));
    z0 = z0/norm(z0,'fro');
  end
  normest = sqrt(sum(Y(:))/numel(Y(:)));
  zinit = normest * z0;
end 

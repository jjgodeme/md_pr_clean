% Implementation of the Mirror Descent (MD) algorithm presented in the paper
% "Phase Retrieval via Mirror Descent: Theory and Algorithms"
% by JJ GODEME & Jalal Fadili
%  Code to reproduce the Phase transition (Monte-Carlo approach)
% The input data are coded diffraction patterns with Real number valued 1D signal
%-n: the size of signal vectors
% -L: is the guess of the relative smoothness constant
% - solving 3rd degre polynom when n become lager (Ne pas oublier augementer la puissance iterer  when n is larger)
% -D: number of masks

clear all;

% Global Parameters.
npower_iter= 50;		% Number of iteratios sin power iteration.
nsample    = 50;		% Number of MC samples.
nlist	   = [8:10:128];	% Signal length.
mlist	   = [128:50:1280];	% Vector of number of measurements.
niter	   = [500 500 500 700];	% Number of iterations of our algorithms.
lambda	   = 1/2;
rho	   = 1/100;	
L	   = 3+rho;   		% Relative smoothness constant (Gaussian model).                                                                                                       
sigma	   = lambda-rho; 	% Local relative strong convexity constant (Gaussian model). 
outerror=zeros(length(nlist),length(mlist),4);	% Output error.
dim=1;
SignalType = 'random';
MeasType   = 'real gaussian';
% Make signal (random_ square, smooth rough or random)
for in =1:length(nlist)
    n=nlist(in);
    x = create_signal(n,SignalType,dim);
    x=x./norm(x,'fro');
    % Function for recording the  errors
    recerror = @(z)min(norm(x-z,'fro'),norm(x+z,'fro'));
    
    for im=1:length(mlist)
        m=mlist(im);
        samplerror=zeros(nsample,4);
        for ta=1:nsample
            % Create operator for the problem
            [A,At] = create_Operator(n,m,dim,MeasType);
            
            %Generate  data measurements
            
            Y = abs(A(x)).^2;% + randn(n,D);
            
            % Objective function, entropy and their derivatives
            
            f = @(z) norm((abs(A(z)).^2-Y),2)^2/(4*m);
            gradf   = @(z) At((abs(A(z)).^2-Y).*A(z))/m;
            psi     = @(z) 1/2*norm(z,'fro')^2+1/4*norm(z,'fro')^4;
            gradpsi = @(z) z.*(norm(z,'fro')^2+1);
            %Spectral initilization for the problem
            
            z0 = randn(n,1);
            [zinit,normest]= Spectral_init(z0,Y,A,At,npower_iter);
            
            % Performing the Mirror descent with backtracking without spectral init.
            options_MDB.niter = niter(4);
            options_MDB.Lbnd  = L;
            options_MDB.kappa = 1e-2;
            options_MDB.eta   = 1.2;
            zMDB=perform_MDBT(randn(n,1),gradf,f,gradpsi,psi,options_MDB);
	    
            %Performing the Mirror descent with spetral init.
            options_MDBI.niter=niter(3);
            options_MDBI.kappa=1e-2;
            options_MDBI.Lbnd=L;
            options_MDBI.eta=1.2;
            zMDBI=perform_MDBT(zinit,gradf,f,gradpsi,psi,options_MDBI);
            
            % Performing the Wirtinger Flow.
            options_WF.niter = niter(1);
            options_WF.normest = normest;
            zWF=perform_WF(zinit,gradf,options_WF);
            
            % Performing the Polyak subgradient method.
            zSG=perform_SG(zinit,gradf,f,niter(2));
            
            %Change global for printing
            samplerror(ta,:)=[recerror(zWF) recerror(zSG) recerror(zMDB) recerror(zMDBI)];
        end
        parfor i=1:4
	    outerror(in,im,i)=length(samplerror(samplerror(:,i)<1e-5))/nsample;
        end
	progressbar(im+(in-1)*length(mlist),length(nlist)*length(mlist));
    end
    %
end
display('well done');
save(['phase_diagram_' MeasType '.mat']);
affichage_file(['phase_diagram_' MeasType '.mat'], 'monte-carlo', MeasType);

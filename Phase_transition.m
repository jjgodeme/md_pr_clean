% Code to reproduce  the Phase transitions (Probability of recovery) of the paper "Phase retrieval via Mirror Descent"
% "Phase Retrieval via Mirror Descent: Theory and Algorithms"
% by JJ GODEME & Jalal Fadili
% The input data are coded diffraction patterns  with Real number valued 1D signal
%-n: the size of signal vectors
% -L: is the guess of the Lipschitz  smad
% - solving 3rd degre polynom when n becomes larger (Ne pas oublier augementer la puissance iterer )
% -D: number of masks

clear;
clearvars; 

% Global Parameters.
npower_iter= 50;		% Number of iteratios in power iteration.
nsample    = 50;		% Number of MC samples.
n	   = 2^3;		% Signal length.
mlist	   = n*[1:14];		% Vector of number of measurements.
niter	   = 500;		% Number of iterations of our algorithms.
lambda	   = 1/2;
rho	   = 1/100;	
MeasType   = 'real gaussian';
L	   = 3+rho;   		% Relative smoothness constant (Gaussian model).                      
MeasType   = 'real cdp';
L	   = 0.8;   		% Relative smoothness constant (CDP model).

sigma	   = lambda-rho; 	% Local relative strong convexity constant (Gaussian model). 
%L=2*(1+rho)^3;			% Relative smoothness constant (CDP model).            
%sigma=lambda-2*rho;		% Local relative strong convexity constant (CDP model).
outerror   = zeros(length(mlist),4);	% Output error.
dim	   = 1;
SignalType = 'random';
x	   = create_signal(n,SignalType,dim);
x	   = x./norm(x,'fro');
% Function for recording the  errors
recerror   = @(z)min(norm(x-z, 'fro'),norm(x+z, 'fro'));

for im=1:length(mlist)
    m=mlist(im);
    samplerror=zeros(nsample,5);
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
        options_MDB.niter = niter;
        options_MDB.Lbnd  = L;
        options_MDB.kappa = 1e-2;
        options_MDB.eta   = 1.2;
        zMDB=perform_MDBT(randn(n,1),gradf,f,gradpsi,psi,options_MDB);
        
        %Performing the Mirror descent  with  spectral init.
        options_MDBI.niter = niter*1.5;
        options_MDBI.kappa = 1e-2;
        options_MDBI.Lbnd  = L;
        options_MDBI.eta   = 1.2;
        zMDBI=perform_MDBT(zinit,gradf,f,gradpsi,psi,options_MDBI);
        
        % Performing the Wirtinger Flow
        options_WF.niter = niter;
        options_WF.normest = normest;
        zWF=perform_WF(zinit,gradf,options_WF);
        
        % Performing the Polyak subgradient method with l1 objective.
	fpol1 = @(z) norm((abs(A(z)).^2-Y),1)/2;
        gradfpol1   = @(z) At(sign(abs(A(z)).^2-Y).*A(z));
        zSG1=perform_SG(zinit,gradfpol1,fpol1,niter*1.5);
	
	% Performing the Polyak subgradient method with squared l^2 objective.
        zSG2=perform_SG(zinit,gradf,f,niter);
        
        %Change global for display.
        samplerror(ta,:)=[recerror(zWF) recerror(zSG1) recerror(zSG2) recerror(zMDB) recerror(zMDBI)];    
    end
    parfor i=1:5
        outerror(im,i)=length(samplerror(samplerror(:,i)<1e-5))/nsample;
    end
    progressbar(im,length(mlist));
end

disp('well done');
save('san_cdpF3.mat');
affichage_file('san_cdpF3.mat', 'probability');

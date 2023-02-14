% Implementation of the Mirror Descent (MD) algorithm presented in the paper
% "Phase Retrieval via Mirror Descent: Theory and Algorithms"
% by JJ GODEME & Jalal Fadili
% The input data are coded diffraction patterns  with Real number valued 1D signal
%-n: the size of signal vectors
% -L: is the guess of the Lipschitz  smad
% - solving 3rd degre polynom when n become lager (Ne pas oublier augementer la puissance iterer )
% -D: number of masks

clear all;

%Global Parameters of the  code

npower_iter= 50;															    % Number of iterate power.
nsample    = 30;															    % Number of sample
nlist      = 2.^[3:8];
mlist      = 1:10:2*floor(nlist(end)*log(nlist(end))); 											% Number of measurements.
niter=700;
delta=1/15;																	    % Number of iterations of our algorithms.
L=3+delta;																	          % Choice smooth apdatable coeficients
outerror=zeros(length(nlist),length(mlist),4);												%  Output error 

% Make signal (random_ square, smooth rough or random)
dim=1;
SignalType = 'smooth rough';

for in =1:length(nlist)
	n=nlist(in);
	x = create_signal(n,SignalType,dim);
	for im =1:length(mlist)
		m=mlist(im);
		samplerror=zeros(nsample,4);
		for ta=1:nsample
		% Create operator for the problem
			[A,At] = create_Operator(n,m,dim,'real gaussian');

		%Generate  data measurements

			Y = abs(A(x)).^2;% + randn(n,D);

		% Objective function, entropy and their derivatives

			f       = @(z) norm((abs(A(z)).^2-Y),2)^2/(4*m);
			gradf   = @(z) At((abs(A(z)).^2-Y).*A(z))/m;
			psi     = @(z) 1/2*norm(z,'fro')^2+1/4*norm(z,'fro')^4;
			gradpsi	= @(z) z.*(norm(z,'fro')^2+1);

		%Spectral initilization for the problem

			z0 = randn(n,1);
			[zinit,normest]= Spectral_init(z0,Y,A,At,npower_iter);

		% Function for recording the  errors

			recerror = @(z)min(norm(x-z, 'fro'),norm(x+z, 'fro'))/norm(x,'fro');

		% Performing the Mirror descent

			options_MD.niter = niter;
			options_MD.kappa = 0;
			options_MD.L     = L;
			zMD=perform_MD(zinit,gradf,gradpsi,options_MD);

		% Performing the Mirror descent with backtracking

			options_MDB.niter = niter;
			options_MDB.Lbnd  =L;
			options_MDB.kappa = 0;
			options_MDB.eta   = 1.2;
			zMDB=perform_MDBT(zinit,gradf,f,gradpsi,psi,options_MDB);

		% Performing the Wirtinger Flow

			options_WF.niter = niter;
			options_WF.normest = normest;
			zWF=perform_WF(zinit,gradf,options_WF);

	 	% Performing the Polyak subgradient methods

			zSG=perform_SG(zinit,gradf,f,niter);

		% Change global for printing
			samplerror(ta,1)=recerror(zWF);
			samplerror(ta,2)=recerror(zSG);
			samplerror(ta,3)=recerror(zMD);
			samplerror(ta,4)=recerror(zMDB);
			progressbar((in-1)*length(mlist)*nsample+(im-1)*nsample+ta,length(nlist)*length(mlist)*nsample)
		end
		for i=1:4
			%outerror(in,im,i)=length(samplerror(samplerror(:,i)<1e-5))/nsample;
			outerror(in,im,i)=mean(samplerror(find(~isnan(samplerror(:,i))),i));
		end
	end
end
 save('Gauss_monte.mat','outerror');
%{
mlist=1:10:2*floor(256*log(256));
        mlist=mlist(1:D);mlist=floor(mlist/128);
        col={'b--+','r--+','k--o','m--o'};
        data=outerror(:,1:D,:);
        for i=1:4
            plot(floor(mlist/128),data(5,:,i),col{i});
            axis;
            hold on;

%Dlist=2*Dlist;
            %Dlist=Dlist(1:floor(length(data)/2));
            %data=data(Dlist);
            %plot(Dlist,data,col{i});
            plot(floor(Dlist/13),data,col{i});


%}
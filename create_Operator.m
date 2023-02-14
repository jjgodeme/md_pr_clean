function [A, At]=create_Operator(n,m,dim,ProblemType)
% This function create the operaor A  for Gaussian and CDP Type
% measurements
if dim==1
    switch lower(ProblemType)
        case 'real gaussian'
            M =1/sqrt(2)*randn(m,n);
            A = @(x) M*x;
            At= @(y) transpose(M)*y;
        case 'real cdp'
            D=(m/n<1)+floor(m/n)*(m/n>=1);
            Masks= reshape(randsample([-1 0 1],n*D,true,[1/4 1/2 1/4]),n,D);
            % Make linear operators;
             A  = @(x)fft(Masks .* repmat(x,[1 D]));	% Input is n x 1 signal, output is n x L array.
             At = @(y)real(sum(Masks .* ifft(y), 2))*n;	% Input is n x L array, output is n x 1 signal.
        otherwise
            disp('Unknown Operator Type')
    end
elseif dim==2
    switch lower(ProblemType)
        case 'real gaussian'
            M=1/sqrt(2)*randn(m,n,n);
            A = @(x) M*x;
            At= @(y) transpose(M)*y;
        case 'real cdp'
            D=(m/n<1)+floor(m/n)*(m/n>=1);
            Masks = zeros(n,n,D);  %Storage for L masks, each of dim n x n
            % Sample phases: each symbol in alphabet {1, 0,-1} has equal prob.
            for ll = 1:D
                Masks(:,:,ll) = reshape(randsample([-1 0 1],n^2,true,[1/4 1/2 1/4]),n,n);
            end
            A = @(I)  fft2(Masks.* reshape(repmat(I,[1 D]), size(I,1), size(I,2), D));
            At = @(Y) real(sum(Masks .* ifft2(Y), 3))*(n^2);
    end
else
     disp('Dimension should be 1 or 2 ');
end

end

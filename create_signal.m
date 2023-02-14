function x= create_signal(n,SignalType,dim)
% This function create an initial signal for the problem of phase
% retrieval mainly take an integer and return a signal depending on your
% choice in  2 D it can be either a  rough surfaceor an image
% dim stand for dimension dim should be equal to 1 or 2

if  dim==1
    switch lower(SignalType)
        case 'random square'
            k = 10;
            ind = randperm(n,k);
            x = zeros(n,1);
            x(ind) = rand(k,1)-1/2;
            x = cumsum(x);
        case 'smooth rough'
            t = linspace(-1/2,1/2,n)';
            s=0.01;
            h=exp(-t.^2/(2*s^2))/sqrt(2*pi*s^2);
            x=real(ifft(fft(rand(n,1)-1/2).*fft(h)))/n;
        otherwise
            x=randn(n,1);
    end
elseif dim==2
    switch lower(SignalType)
        case 'smooth rough'
            t = linspace(-1/2,1/2-1/n,n)';
            [T1,T2] = meshgrid(t,t);
            s=0.01;
            h=exp(-(T1.^2+T2.^2)/(2*s^2))/(2*pi*s^2);
            x=10*real(ifft2(fft2(rand(n,n)-1/2).*fft2(h)))/n^2;
        case 'phenix'
            x=double(imread('phenix.jpg'));
            x=imresize(mean(x,3),[n n]);
        otherwise
            x=randn(n,n);
    end
else
    disp('Error: Dimension should be 1 or 2');
end

end

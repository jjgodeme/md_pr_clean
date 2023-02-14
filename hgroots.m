function z=hgroots(P)
P21=P(2)/P(1);
P31=P(3)/P(1);
P41=P(4)/P(1);
p=-P21^2/3+P31;
q=(2*P21^3)/27-(P21*P31)/3+P41;
delta= q^2+4/27*p^3; 
if delta>0
 u=nthroot((-q+sqrt(delta))/2,3);
 if -q-sqrt(delta)/2<0
     v=-nthroot((q+sqrt(delta))/2,3);
 else     
     v=nthroot((-q-sqrt(delta))/2,3);
 end 
 z=u+v-P21/3;
elseif delta <0
    u=nthroot((-q+1i*(sqrt(-2*delta)))/2,3);
    z=real(u+conj(u)-P21/3);
else 
    z=(3*q/p-P21/3);
end 
end 

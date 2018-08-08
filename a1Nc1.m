function [h,Time]=a1Nc1(eps,N)%

T = 10^6;
 
delta=1/(2*(N+1)); %%Fix the delta such that delta < 1/N

for j=1:N+1
    h(1,j)=input('enter position of layers in ascending order:   ');
end
h0(1,1)=-h(1,1);
h(1,N+2)=2-h(1,N+1);
%%define Mj's
m(1,1)=(h0(1,1)+h(1,1))/2;
for j=2:N+2
    m(1,j)=(h(1,j-1)+h(1,j))/2;
end
s=eps/10;%%Step Size
%%Define hyperbolic function( for f(u)=u^3-u)
A=sqrt(2);
K=4;

 %%Iteration loop
 dt=input('Enter time step:  ');
 
 for k=1:T;
     for j=1:N+1
    if h(k,j)<=0 
        h(k,j)=0;
    end
     if h(k,j)>=1
            h(k,j)=1;
     end
     
     end
       rk1(k,1)=c(h(k,1) , h(k,2) , h(k,3) , eps)*dt;
       rk1(k,2)=r(h(k,1),h(k,2),h(k,3),h(k,4),eps)*dt;
        for j=3:N-1 
       rk1(k,j)=w(h(k,j-2),h(k,j-1),h(k,j),h( k,j+1),h(k,j+2),eps)*dt;
        end
        rk1(k,N)=q(h(k,N-2),h(k,N-1),h(k,N),h( k,N+1),eps)*dt;
        rk1(k,N+1)=o(h(k,N-1),h(k,N),h(k,N+1),eps)*dt;
        
        
       rk2(k,1)=c(h(k,1)+rk1(k,1)/2,h(k,2)+rk1(k,2)/2,h(k,3)+rk1(k,3)/2,eps)*dt;
       rk2(k,2)=r(h(k,1)+rk1(k,1)/2,h(k,2)+rk1(k,2)/2,h(k,3)+rk1(k,3)/2,h(k,4)+rk1(k,4)/2,eps)*dt;
       for j=3:N-1
            rk2(k,j)=w(h(k,j-2)+rk1(k,j-2)/2,h(k,j-1)+rk1(k,j-1)/2,h(k,j)+rk1(k,j)/2,h(k,j+1)+rk1(k,j+1)/2,h(k,j+2)+rk1(k,j+2)/2,eps)*dt;
       end
       rk2(k,N)=q(h(k,N-2)+rk1(k,N-2)/2,h(k,N-1)+rk1(k,N-1)/2,h(k,N)+rk1(k,N)/2,h(k,N+1)+rk1(k,N+1)/2,eps)*dt;
        rk2(k,N+1)=dt*o(h(k,N-1)+rk1(k,N-1)/2,h(k,N)+rk1(k,N)/2,h(k,N+1)+rk1(k,N+1)/2,eps);
        
        %third iteration
        
       rk3(k,1)=c(h(k,1)+rk2(k,1)/2,h(k,2)+rk2(k,2)/2,h(k,3)+rk2(k,3)/2,eps)*dt;
       rk3(k,2)=r(h(k,1)+rk2(k,1)/2,h(k,2)+rk2(k,2)/2,h(k,3)+rk2(k,3)/2,h(k,4)+rk2(k,4)/2,eps)*dt;
       for j=3:N-1
       rk3(k,j)=w(h(k,j-2)+rk2(k,j-2)/2,h(k,j-1)+rk2(k,j-1)/2,h(k,j)+rk2(k,j)/2,h(k,j+1)+rk2(k,j+1)/2,h(k,j+2)+rk2(k,j+2)/2,eps)*dt;
       end
        rk3(k,N)=q(h(k,N-2)+rk2(k,N-2)/2,h(k,N-1)+rk2(k,N-1)/2,h(k,N)+rk2(k,N)/2,h(k,N+1)+rk2(k,N+1)/2,eps)*dt;
       rk3(k,N+1)=dt*o(h(k,N-1)+rk2(k,N-1)/2,h(k,N)+rk2(k,N)/2,h(k,N+1)+rk2(k,N+1)/2,eps);
       
       %4th iteration
       
       rk4(k,1)=c(h(k,1)+rk3(k,1),h(k,2)+rk3(k,2),h(k,3)+rk3(k,3),eps)*dt;
       rk4(k,2)=r(h(k,1)+rk3(k,1),h(k,2)+rk3(k,2),h(k,3)+rk3(k,3),h(k,4)+rk3(k,4),eps)*dt;
       for j=3:N-1
       rk4(k,j)=w(h(k,j-2)+rk3(k,j-2),h(k,j-1)+rk3(k,j-1),h(k,j)+rk3(k,j),h(k,j+1)+rk3(k,j+1),h(k,j+2)+rk3(k,j+2),eps)*dt;
       end
       rk4(k,N)=q(h(k,N-2)+rk3(k,N-2),h(k,N-1)+rk3(k,N-1),h(k,N)+rk3(k,N),h(k,N+1)+rk3(k,N+1),eps)*dt;
       rk4(k,N+1)=dt*o(h(k,N-1)+rk3(k,N-1),h(k,N)+rk3(k,N),h(k,N+1)+rk3(k,N+1),eps);
      % calculating hj
       h(k+1,1)=h(k,1)+(1/6)*(rk1(k,1)+2*rk2(k,1)+2*rk3(k,1)+rk4(k,1));
        h(k+1,2)=h(k,2)+(1/6)*(rk1(k,2)+2*rk2(k,2)+2*rk3(k,2)+rk4(k,2));
     for j=3:N
        h(k+1,j)=h(k,j)+(1/6)*(rk1(k,j)+2*rk2(k,j)+2*rk3(k,j)+rk4(k,j));
     end
      h(k+1,N)=h(k,N)+(1/6)*(rk1(k,N)+2*rk2(k,N)+2*rk3(k,N)+rk4(k,N));
     h(k+1,N+1)=h(k,N+1)+(1/6)*(rk1(k,N+1)+2*rk2(k,N+1)+2*rk3(k,N+1)+rk4(k,N+1));
     a=k;
     Time= a*dt;
     if h(k+1,1)+eps>=h(k+1,2);
    
         break
     elseif h(k+1,2)+eps>=h(k+1,3);
         
         break
         elseif h(k+1,3)+eps>=h(k+1,4);
         
         break
       elseif h(k+1,1)<=eps;
         
         break
     elseif h(k+1,4)+eps>=1;
         
         break
     end
         
      
 end

s=eps/10;%%Step Size
%%Define hyperbolic function( for f(u)=u^3-u) 
%Plotting function
%title(['Time is ',num2str(Time),' s'])
subplot(2,2,1)
plot(h(:,1),'k');
title(['h_' num2str(1),' Total Time is ',num2str(Time),' s' ])
xlabel('Time step')
ylabel(['h_' num2str(1)])
 
subplot(2,2,2)
plot(h(:,2),'k');
title(['h_' num2str(2)])
xlabel('Time step')
ylabel(['h_' num2str(2)])

subplot(2,2,3)
plot(h(:,3),'k');
title(['h_' num2str(3)])
xlabel('Time step')
ylabel(['h_' num2str(3)])

subplot(2,2,4)
plot(h(:,4),'k');
title(['h_' num2str(4)])
xlabel('Time step')
ylabel(['h_' num2str(4)])
 
end
 
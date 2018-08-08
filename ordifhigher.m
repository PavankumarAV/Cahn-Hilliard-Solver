function [u,h,Time]=ordifhigher(Num,eps,a,T,Time,b)
%a=Num %b=eps, c=h(1)d=h(2) dt=e , f=T, g=time
h=a; % input from switch% difference between layers
h0=-h(1);% ghost layer left
h(Num+2)=2-h(Num+1); %ghost layer right
m(1)=0; % always constant for [0,1] case
dt=0.01;% initial guess as small as possible
for k=1:T;
    d=abs(r(h(1),h(2),h(3),h(4),eps));
    
    if b==1% general variable for time adaptivity
    
    dt=1/(10000*d);
    elseif b==2;
        dt=Settime(d);
    end% continuous time adaptivity, if movement is too fast, change 2000 to 5000 or 10,000 to see slow movement
    % runge kutta fourth order iteration   
    rk1(1)=c(h(1) , h(2) , h(3) , eps)*dt;
       rk1(2)=r(h(1),h(2),h(3),h(4),eps)*dt;
        for j=3:Num-1 
       rk1(j)=w(h(j-2),h(j-1),h(j),h(j+1),h(j+2),eps)*dt;
        end
        rk1(Num)=q(h(Num-2),h(Num-1),h(Num),h(Num+1),eps)*dt;
        rk1(Num+1)=o(h(Num-1),h(Num),h(Num+1),eps)*dt;
        
        
       rk2(1)=c(h(1)+rk1(1)/2,h(2)+rk1(2)/2,h(3)+rk1(3)/2,eps)*dt;
       rk2(2)=r(h(1)+rk1(1)/2,h(2)+rk1(2)/2,h(3)+rk1(3)/2,h(4)+rk1(4)/2,eps)*dt;
       for j=3:Num-1
            rk2(j)=w(h(j-2)+rk1(j-2)/2,h(j-1)+rk1(j-1)/2,h(j)+rk1(j)/2,h(j+1)+rk1(j+1)/2,h(j+2)+rk1(j+2)/2,eps)*dt;
       end
       rk2(Num)=q(h(Num-2)+rk1(Num-2)/2,h(Num-1)+rk1(Num-1)/2,h(Num)+rk1(Num)/2,h(Num+1)+rk1(Num+1)/2,eps)*dt;
        rk2(Num+1)=dt*o(h(Num-1)+rk1(Num-1)/2,h(Num)+rk1(Num)/2,h(Num+1)+rk1(Num+1)/2,eps);
        
        %third iteration
        
       rk3(1)=c(h(1)+rk2(1)/2,h(2)+rk2(2)/2,h(3)+rk2(3)/2,eps)*dt;
       rk3(2)=r(h(1)+rk2(1)/2,h(2)+rk2(2)/2,h(3)+rk2(3)/2,h(4)+rk2(4)/2,eps)*dt;
       for j=3:Num-1
       rk3(j)=w(h(j-2)+rk2(j-2)/2,h(j-1)+rk2(j-1)/2,h(j)+rk2(j)/2,h(j+1)+rk2(j+1)/2,h(j+2)+rk2(j+2)/2,eps)*dt;
       end
        rk3(Num)=q(h(Num-2)+rk2(Num-2)/2,h(Num-1)+rk2(Num-1)/2,h(Num)+rk2(Num)/2,h(Num+1)+rk2(Num+1)/2,eps)*dt;
       rk3(Num+1)=dt*o(h(Num-1)+rk2(Num-1)/2,h(Num)+rk2(Num)/2,h(Num+1)+rk2(Num+1)/2,eps);
       
       %4th iteration
       
       rk4(1)=c(h(1)+rk3(1),h(2)+rk3(2),h(3)+rk3(3),eps)*dt;
       rk4(2)=r(h(1)+rk3(1),h(2)+rk3(2),h(3)+rk3(3),h(4)+rk3(4),eps)*dt;
       for j=3:Num-1
       rk4(j)=w(h(j-2)+rk3(j-2),h(j-1)+rk3(j-1),h(j)+rk3(j),h(j+1)+rk3(j+1),h(j+2)+rk3(j+2),eps)*dt;
       end
       rk4(Num)=q(h(Num-2)+rk3(Num-2),h(Num-1)+rk3(Num-1),h(Num)+rk3(Num),h(Num+1)+rk3(Num+1),eps)*dt;
       rk4(Num+1)=dt*o(h(Num-1)+rk3(Num-1),h(Num)+rk3(Num),h(Num+1)+rk3(Num+1),eps);
      % calculating hj
       hn(1)=h(1)+(1/6)*(rk1(1)+2*rk2(1)+2*rk3(1)+rk4(1));
        hn(2)=h(2)+(1/6)*(rk1(2)+2*rk2(2)+2*rk3(2)+rk4(2));
     for j=3:Num
        hn(j)=h(j)+(1/6)*(rk1(j)+2*rk2(j)+2*rk3(j)+rk4(j));
     end
      hn(Num)=h(Num)+(1/6)*(rk1(Num)+2*rk2(Num)+2*rk3(Num)+rk4(Num));
     hn(Num+1)=h(Num+1)+(1/6)*(rk1(Num+1)+2*rk2(Num+1)+2*rk3(Num+1)+rk4(Num+1));
        hn(Num+2)=2-hn(Num+1);
        s=eps/10;%%Step Size
        m(1)=0;
   for j=2:Num+2
             m(j)=(h(j-1)+h(j))/2;
   end
   x1=0:s:m(2);
   x2=m(2):s:m(3);
   x3=m(3):s:m(4);
   x4=m(4):s:1;
   u1=((-1)^(2)).*tanh((x1-hn(1))/eps);
   u2=((-1)^(2+1)).*tanh((x2-hn(2))/eps);
   u3=((-1)^(3+1)).*tanh((x3-hn(3))/eps);
   u4=((-1)^(4+1)).*tanh((x4-hn(4))/eps);
Time=Time+dt;
xp=[x1, x2, x3, x4];
up=[u1,u2,u3,u4];
plot(xp,up);
xlabel('x------>');
        ylabel('u');
title([' Total Time is ',num2str(Time),' s' ])
drawnow;
h=hn;
u=up;
if h(1)>=h(2)-3*eps;
    break
elseif h(2)>=h(3)-3*eps; % making dependent on eps
     break
elseif h(3)>=h(4)-3*eps;
    break
elseif h(4)>=1-3*eps;
    break
end
       
 end 
end % returns the value Time and h=[h(1), h(2)]

function [h,Time]=a1c2e2(N,eps,e)
Num=N; % define the number of layers(Num+1 layers)
 for j=1:N+1
    h(j)=input('enter position of layers in ascending order:   ');
  end
T=10000;
hdif=h(2)-h(1);% difference between layers
h0=-h(1);% ghost layer left
h(3)=2-h(2); %ghost layer right
m(1)=0; 
if e==1% always constant for [0,1] case
delt=input('Enter time step:   ');
end% initial guess as small as possible
Time=0;
for k=1:T;
    d=abs(f(h(1),h(2),eps)); 
    if e==1
        dt=delt; %constant time
        %dt=Settime(d); % you can also use this for discountinuous
        %adaptivity
    elseif e==2
        dt=1/(10000*d); % continuous time adaptivity, if movement is too fast, change 20000 to 50000 to see slow movement
    end
    
     
    % runge kutta fourth order iteration   
    rk1(1)=f(h(1),h(2),eps)*dt; 
       rk2(1)=dt*f(h(1)+rk1(1)/2,h(2)+rk1(1)/2,eps);
       rk3(1)=dt*f(h(1)+rk2(1)/2,h(2)+rk2(1)/2,eps);
       rk4(1)=dt*f(h(1)+rk3(1),h(2)+rk3(1),eps);
       hnext(1)=h(1)+(1/6)*(rk1(1)+2*rk2(1)+2*rk3(1)+rk4(1));
       hnext(2)=h(1)+hdif;
       h(1)=hnext(1); % new position oh h(1)
       h(2)=hnext(2); % new position of h(2)
      h0(1)=-h(1); % new positin of Ghost left
      h(3)=2-h(2);% new positin of Ghost right
      Time=Time+dt;  % count time 
      m(2)=(h(1)+h(2))/2; % new position
      m(3)=1;
        
        s=eps/10;%%Step Size
        
        %%Define hyperbolic function( for f(u)=u^3-u) 
         x1=m(1):s:m(2);
         x2=m(2):s:1;
        u1=((-1)^(1+1)).*tanh((x1-h(1))/eps);
        u2=((-1)^(2+1)).*tanh((x2-h(2))/eps);
        u=[u1 u2];
        x=[x1 x2];
        plot(x,u)
        xlabel('x------>')
       ylabel('u')
        axis([0 1 -1 1])
        title(['Time is ',num2str(Time),' s'])
        drawnow
        %loop break condition, if the layers are expected go towards 1,
        %then change the condition to h(2)=>1-2*eps
        if h(1)<=2*eps % making dependent on eps
            break
        elseif h(2)>=1-2*eps % making dependent on eps
            break
        end
       
 end 
end % returns the value Time and h=[h(1), h(2)]

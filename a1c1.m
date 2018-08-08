function [h]=a1c1(eps,N)
T=input('Enter the time scale:   ');
h=zeros(T,N+2);
  L=zeros(T,N+2);
  alpha=zeros(T,N+2);
  m=zeros(T,N+2);
  cons=zeros(T,N+2);
  g=zeros(T,N+2);
delta=1/(2*(N+1));   

%%input the Integer such that there are N+1 transition layers in the solution.
%% for example if there are 2 Transition layer, then enter 1.
for j=1:N+1
    h(1,j)=input('enter position of layers in ascending order:   ');
end
h0(1,1)=-h(1,1);
h(1,N+2)=2-h(1,N+1);
m(1,1)=(h0(1,1)+h(1,1))/2;
for j=2:N+2
    m(1,j)=(h(1,j-1)+h(1,j))/2;
end
s=eps/10;%%Step Size
%%Define hyperbolic function( for f(u)=u^3-u)
A=sqrt(2);
K=4;
L(1,1)=h(1,1)-h0(1,1);
for j=2:N+2
    L(1,j)=h(1,j)-h(1,j-1);
end
 %initialize aplha values and constants
 c1=-A/eps;
 for j=1:N+2
     alpha(1,j)=exp(c1*L(1,j));
 end
 cons(1,1)=(16/(L(1,2)))*(alpha(1,3)-alpha(1,1));
 cons(1,N+1)=(16/(L(1,2)))*(alpha(1,3)-alpha(1,1));
 for j=2:N
     cons(1,j)=(16/(L(1,j)))*(alpha(1,j+1)-alpha(1,j-1))+(16/(L(1,j+1)))*(alpha(1,j+2)-alpha(1,j)) ;
 
 end
 %%Iteration loop 
 for k=1:T;
    if h(k,1)<=0 
        h(k,1)=0;
        break
    end
        
       g(k,1)=cons(k,1)*0.1;
       g(k,2)=cons(k,2)*0.1;
      h(k+1,1)=g(k,1)+ h(k,1);
       h(k+1,2)=g(k,2)+ h(k,2);
      h0(k+1,1)=-h(k+1,1);
      h(k+1,3)=2-h(k+1,2);
      L(k+1,1)=2*h(k+1,1) ;
      L(k+1,2)=h(k+1,2)-h(k+1,1);
      L(k+1,3)=h(k+1,3)-h(k+1,2);
      alpha(k+1,1)=exp(c1*L(k+1,1));
      alpha(k+1,3)=exp(c1*L(k+1,3));
      cons(k+1,1)=(16/(L(k+1,2)))*(alpha(k+1,3)-alpha(k+1,1));
      cons(k+1,2)=(16/(L(k+1,2)))*(alpha(k+1,3)-alpha(k+1,1));
      
 end
  
s=eps/10;%%Step Size
%%Define hyperbolic function( for f(u)=u^3-u) 
for k=1:T
for j=1:N+1
    x{k,j}=m(k,j):s:m(k,j+1); % define Array to store values of x as cut-off function
 u{k,j}=((-1)^(j+1)).*tanh((x{k,j}-h(k,j))/eps); %Piercing solutions in array and stored in u  %end
end
end
  

 plot(h(:,1))
  hold on
  plot(h(:,2))
   xlabel('t------>')
       ylabel('h1, h2')
  hold off

 
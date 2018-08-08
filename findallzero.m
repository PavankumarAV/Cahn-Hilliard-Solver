function [N,Zn]=findallzero(x,y,m)
% x=xp;
% y=u;
% m=length(u);
 N=0;
for i=1:m-1
    if y(i)*y(i+1)<=0
        z(i)=x(i)-((x(i+1)-x(i))/(y(i+1)-y(i)))*y(i); 
        N=N+1;
    end
end
Zn=z(z~=0); 

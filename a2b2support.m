function [u]=a2b2support(N,h,eps)

delta=1/(2*(N+1)); %%Fix the delta such that delta < 1/N

h0=-h(1);
h(N+2)=2-h(N+1);
%%define Mj's
m(1)=0;

for j=2:N+2
    m(j)=(h(j-1)+h(j))/2;
end
s=eps/10;%%Step Size 
for j=1:N+1
    x{j}=m(j):s:m(j+1); 
 up{j}=((-1)^(j+1)).*tanh((x{j}-h(j))/eps);
end

u=up{1};
xp=x{1};
for j=2:N+1
    u=[u , up{j}];
    xp=[xp, x{j}];
end
x=xp;
end
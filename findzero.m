function z=findzero(x,y,m)

for i=1:m-1
    if y(i)*y(i+1)<=0
        z=x(i)-((x(i+1)-x(i))/(y(i+1)-y(i)))*y(i); 
    end
end


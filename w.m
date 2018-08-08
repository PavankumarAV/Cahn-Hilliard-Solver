function w=w(x,y,z,s,g,eps)
w=(1/4*(z-y))*(exp(-sqrt(2)*((s-z))/eps)-exp(-sqrt(2)*((y-x))/eps))+(1/4*(s-z))*(exp(-sqrt(2)*((g-s))/eps)-exp(-sqrt(2)*((z-y))/eps));
end 

function q=q(x,y,z,s,eps)
q=(1/4*(z-y))*(exp(-sqrt(2)*((s-z))/eps)-exp(-sqrt(2)*((y-x))/eps))+(1/4*(s-z))*(exp(-sqrt(2)*((2-2*s))/eps)-exp(-sqrt(2)*((z-y))/eps));
end 

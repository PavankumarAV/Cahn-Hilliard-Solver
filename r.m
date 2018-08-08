function r=r(x,y,z,g,eps)
r=(1/4*(y-x))*(exp(-sqrt(2)*(z-y)/eps)-exp(-sqrt(2)*(2*x)/eps))+(1/4*(z-y))*(exp(-sqrt(2)*((g-z))/eps)-exp(-sqrt(2)*((y-x))/eps));
end 
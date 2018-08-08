function c=c(x,y,z,eps)
c=(1/4*(y-x))*(exp(-sqrt(2)*((z-y))/eps)-exp(-sqrt(2)*((2*x))/eps));
end 
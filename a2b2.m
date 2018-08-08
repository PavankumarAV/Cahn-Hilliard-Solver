function [Time]=a2b2(eps,b)
if b==1
     h = eps/10;
        M=1/h;
         delx = 1/(M-1);
         x = (0:delx:1)';
         T=100000;
         u=input('Press initial condition for example cos(5*x), sin(50*x):  ')';
U = u';
elseif b==3
    C=input('Enter the constant mean value Between -1 to +1:   ');
        h = eps/10;
        M=1/h;
         delx = 1/(M-1);
         x = (0:delx:1)';
        U=C+ 0.001*rand(M,1);
elseif b==2
    N=input('Enter the number of N+1 transition layers:   ');
   for j=1:N+1
   h(j)=input('enter position of layers in ascending order:  ');
   end
    [u]=a2b2support(N,h,eps);
    M=length(u);
        delx = 1/(M-1);
        delx2 = delx^2;
        x = (0:delx:1)';
    U=u';
end
dt= input('Enter positive time step, Less than 0.001:  ');
delx2 = delx^2;
delt = dt;
ntmax = 10^6;
epsilon = eps;
eps2 = epsilon^2;
a = 2;
lam1 = delt/delx2;
Time=0;
lam2 = eps2*lam1/delx2;
Leig  = (((2*cos(pi*(0:M-1)'/(M-1)))-2));
 % scaled eigenvalues of stabilized CH update matrix
CHeig = ones(M,1) - (a*lam1*Leig) + (lam2*Leig.*Leig);
% scaled eigenvalues of the laplacian
Seig = lam1*Leig;
% random initial conditions
 hat_U = dct2(U);
% main loop
        
        for it=1:ntmax

          Time = Time+delt;
        % compute the shifted nonlinear term
          fU = (U.*U.*U) - ((1+a)*U);
        % compute the right hand side in tranform space
          hat_rhs = hat_U + (Seig.*dct2(fU));
        % compute the updated solution in tranform space
          hat_U = hat_rhs./CHeig;
        % invert the cosine transform
          U = idct2(hat_U);
          plot(x,U);
          xlabel('x------>')
       ylabel('u')
          title([' Total Time is ',num2str(Time),' s' ])
          axis([0 1 -1 1])
          drawnow;% play the movie
                
        end
end

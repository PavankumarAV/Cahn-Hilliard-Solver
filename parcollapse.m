function [u,h,Time]=parcollapse(Num,eps,a,T,Time,z)
           dt=0.0001;
           h=a;
        %a=Num %b=eps, c=h(1)d=h(2) dt=e , f=T, g=tim
       u=z;
        N=length(u); %for mesh
        delx = 1/(N-1); %dx
        delx2 = delx^2;
        x = (0:delx:1)'; % grid generation
        delt = dt; % dt
        ntmax = 1200;% max number of iteration to run
        epsilon = eps; %eps
        eps2 = epsilon^2;
        a = 2; % pre-conditioner value
        lam1 = delt/delx2; 
        lam2 = eps2*lam1/delx2;
        Leig  = (((2*cos(pi*(0:N-1)'/(N-1)))-2));
        % scaled eigenvalues of stabilized CH update matrix
        CHeig = ones(N,1) - (a*lam1*Leig) + (lam2*Leig.*Leig);
        % scaled eigenvalues of the laplacian
        Seig = lam1*Leig;
        % random initial conditions
        U =  u';
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
          plot(x,U)
          xlabel('x------>');
        ylabel('u');
          title([' Total Time is ',num2str(Time),' s' ])
          axis([0 1 -1 1])
          drawnow;% play the movie
          u=U;
          [N,Zn]=findallzero(x,u,length(u));
          h=Zn;
        end 
end


function [h,Time]=pardif(Num,eps,a,b,T,Time)         
h(1)=a;
            h(2)=b;
            hdif=h(2)-h(1);
        %a=Num %b=eps, c=h(1)d=h(2) dt=e , f=T, g=tim
        h0=-h(1);
        h(3)=2-h(2);
        %%define Mj's
        m(1)=0;
        for j=2:Num+2
            m(j)=(h(j-1)+h(j))/2;
        end
        s=eps/10;%%Step Size
        
        %%Define hyperbolic function( for f(u)=u^3-u) 
        for j=1:Num+1
             x{j}=m(j):s:m(j+1); % define Array to store values of x as cut-off function
         u{j}=((-1)^(j+1)).*tanh((x{j}-h(j))/eps);
        end

        %Plotting function
        %     for i=1:Num+1
        %        plot(x{i},u{i},'b')
        %        xlabel('x------>')
        %        ylabel('u^h(x)')
        %     hold on
        % end
         u=[u{1} u{2}];



        N=length(u);
        delx = 1/(N-1);
        delx2 = delx^2;
        x = (0:delx:1)';
        delt = dt;
        ntmax = T;
        epsilon = eps;
        eps2 = epsilon^2;
        a = 2;
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
                h(1)=findzero(x,U, length(u));
                h(2)=h(1)+hdif;
                if Time>=0.10
                    break
                end
        end
end


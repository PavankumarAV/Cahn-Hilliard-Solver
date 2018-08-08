clc
close all;
clear all;
eps= input('Enter epsilon value:    ');
a= input('For working with ODE--> press 1; For PDE--> press 2; for switching scheme--> press 3 :   ');
if a==1
   c=input('Do you wish to see graph of movement of layer position--> press 1; for movie--> press 2:   ');
    N=input('Enter the number of N-1 transition layers:    ');
    if N==1
        if c==1
       [h]=a1c1(eps,N);% rkmethod 3 layers  
        elseif c==2
           e=input('For constant time step--> Press 1; For continuous adaptivity--> Press 2:   ');
               [h,Time]=a1c2e2(N,eps,e);%continuous adaptivity
         end
   elseif N>2
    if c==1
       [h,Time]=a1Nc1(eps,N) ;% modifiesrk.m  
    elseif c==2
        e=input('For constant time step--> Press 1; For continuous adaptivity--> Press 2:   ');
               [h, Time]=a1Nc2e2(N,eps,e);%continuous adaptivity 4 layer

    end
    end
elseif a==2
    b=input('Initial condition:for Typing user defined condition press 1; for Layered structure Press 2; for Random initial condition around some constant c press 3:   ');
    [Time]=a2b2(eps,b);
elseif a==3
    N=input('Enter the number of N+1 transition layers: N=1 and N=3 are available:    ');
    for j=1:N+1
    h(j)=input('enter position of layers in ascending order:   ');
    end
     b=input('For continuous adaptivity press 1, for discontinuous adaptivity press 2:   ');
   [T1, Time]=a3N3(eps, N, h, b); %a3n1
     
end
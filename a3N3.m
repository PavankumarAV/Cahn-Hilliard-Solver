function [T1, Time]=a3N3(eps, N, h, b)
format long %saves 16 digits after point
Num = N; %Number of Transition layers
T =10000;% Number of Steps at each level % takes the input from User for epsilon value
A=sqrt(2);% Constant Value
K=4; %Constant Value
dt=0.001;%Initial Time step
Time=0; % End Time( Initially Zero)
if N==3
      [u,h,Time]=ordifhigher(Num,eps,h, T,Time,b); % operation in ODE
   [u,h,Time]= parcollapse(Num,eps,h, T,Time,u); %Switch during Collapse
   Num=Num-2; %collapse

T1=Time; % Save intermediate collpase time
% 2 layers starts
[u,h,Time]=ordif(Num,eps,h,T,Time,b);% System perform orperation in ODE..
[Time]= pardif2(Num,eps,h(1),h(2),T,Time,u);
% during the collapse, system switch into PDE
elseif N==1
    Time=0;
    [u,h,Time]=ordif(Num,eps,h,T,Time,b);% System perform orperation in ODE..
[Time]= pardif2(Num,eps,h(1),h(2),T,Time,u);
end
    

end
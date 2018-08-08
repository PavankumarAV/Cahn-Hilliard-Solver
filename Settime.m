function dt=Settime(d)
       if d<10^-50
           dt=1/(10000*d);
         elseif d<10^-28
           dt=10^25;
           elseif d<10^-27
           dt=10^24;
            elseif d<10^-26
           dt=10^23;
           elseif d<10^-25
           dt=10^22;
           elseif d<10^-24
           dt=10^21;
           elseif d<10^-23
           dt=10^20;
           elseif d<10^-22
           dt=10^19;
           elseif d<10^-21
           dt=10^18;
           elseif d<10^-20
           dt=10^17;
           elseif d<10^-19
           dt=10^16;
           elseif d<10^-18
           dt=10^15;
       elseif d<10^-17
           dt=10^14;
       elseif d<10^-16
           dt=10^13;
       elseif d<10^-15
           dt=10^12;
           elseif d<10^-14
           dt=10^11;
       elseif d<10^-13
           dt=10^10;
      elseif d<10^-12
           dt=10^9;
       elseif d<10^-10
           dt=10^7;
       elseif d<10^-8
           dt=10^5;
       elseif d<10^-6
           dt=1000;
       elseif d<10^-5
           dt=100;
           elseif d<10^-4
           dt=10;
           elseif d<10^-3
           dt=1;
           elseif d<10^-2
           dt=0.01;
       elseif d<10^-1
           dt=0.001;
       end
end
nx=51;
dx=(1/(nx-1));
dt=0.1*dx*dx;
c=dt/(2*dx*dx);
Bi=1;
Bi2=10;
temp_first=zeros(nx,nx,nx);%intial condition 
temp=zeros(nx,nx,nx); 
temp_initial=zeros(nx,nx,nx);
temp_second=zeros(nx,nx,nx);
x=linspace(0,1,nx);%x grid
y=linspace(0,1,nx); %y grid
z=linspace(0,1,nx); %z grid
w=0.5; %relaxation factor
error=1;
error_req=1e-8;
itration_max=200000;  %no of iterations
for time=1:200000
    time=time+1;
    rms=1;
 for it=1:itration_max
        %Gauss seidel method
        rms=0;
        res=0;
        %interior points
        for i=2:nx-1
            for j=2:nx-1
                for k=2:nx-1
                 rhs=(dt/c)+((1-6*c)/c)*temp_first(i,j,k)+temp_first(i-1,j,k)+temp_first(i,j-1,k)+temp_first(i,j,k-1)+temp_first(i+1,j,k)+temp_first(i,j+1,k)+temp_first(i,j,k+1);   
                 lhs=((1+6*c)/c)*temp(i,j,k)-(temp(i-1,j,k)+temp(i,j-1,k)+temp(i,j,k-1)+temp(i+1,j,k)+temp(i,j+1,k)+temp(i,j,k+1));   
                 res=rhs-lhs;
                temp(i,j,k)=temp(i,j,k)+w*res/((1+6*c)/c);
                rms=rms+res*res;
              end
            end
        end
        %exterior points (boundary conditions)
        %bc1
        i=1;
        for j=1:nx
           for k=1:nx
           res1=3*temp(i,j,k)-4*temp(i+1,j,k)+temp(i+2,j,k);
           temp(i,j,k)=temp(i,j,k)+w*res1/(-3);
           rms=rms+res1*res1;
           end
        end
        %bc2
        i=nx;
        for j=2:nx
           for k=2:nx
           res2=-(3+2*dx*Bi)*temp(i,j,k)+4*temp(i-1,j,k)-temp(i-2,j,k);
             temp(i,j,k)=temp(i,j,k)+w*res2/(3+2*dx*Bi);
            rms=rms+res2*res2;
           end
        end
        %bc3
        j=1;
        for i=2:nx
          for k=1:nx
           res3=3*temp(i,j,k)-4*temp(i,j+1,k)+temp(i,j+2,k);
           temp(i,j,k)=temp(i,j,k)+w*res3/(-3);
           rms=rms+res3*res3;
          end
        end
        %bc4
        j=nx;
       for i=2:nx-1
         for k=2:nx
         res4=-(3+2*dx*Bi2)*temp(i,j,k)+4*temp(i,j-1,k)-temp(i,j-2,k);
         temp(i,j,k)=temp(i,j,k)+w*res4/(3+2*dx*Bi2);
          rms=rms+res4*res4;
         end
       end
        %bc5
        k=1;
       for i=2:nx
         for j=2:nx
           res5=3*temp(i,j,k)-4*temp(i,j,k+1)+temp(i,j,k+2);
           temp(i,j,k)=temp(i,j,k)+w*res5/(-3);
           rms=rms+res5*res5;
         end
       end
        %bc6
        k=nx;
        for i=2:nx-1
         for   j=2:nx-1
         res6=-(3+2*dx*Bi)*temp(i,j,k)+4*temp(i,j,k-1)-temp(i,j,k-2);
         temp(i,j,k)=temp(i,j,k)+w*res6/(3+2*dx*Bi);
          rms=rms+res6*res6;
         end
        end
        rms=sqrt(rms/(nx*nx*nx)); %root mean square
        
      if rms<error_req
          break 
      end  
     
 end 
    %end of gauss seidel
    
    error=0;
    for i=1:nx
        for j=1:nx
            for k=1:nx
            error=error+((temp(i,j,k)-temp_first(i,j,k))^2);
           end
        end
    end
            
    error=sqrt(error);
    if error<error_req
       break
    end
   for i=1:nx
        for j=1:nx
            for  k=1:nx
        temp_first(i,j,k)=temp(i,j,k);
         end
        end
   end
    if time == 4000
        for i=1:nx
            for j=1:nx
                for k=1:nx
                    temp_initial(i,j,k)=temp(i,j,k);       
                end
            end
        end
        
    end

    if time == 70000
        for i=1:nx
            for j=1:nx
                for k=1:nx
                    temp_second(i,j,k)=temp(i,j,k);
                end
            end
        end
    
    end
end 
%contour plotting
T1=zeros(nx,nx);
T2=zeros(nx,nx);
T3=zeros(nx,nx);
[X,Y] = meshgrid(x,y);

% plots at z=0
for i=1:nx
    for j=1:nx
        T1(i,j) = temp_initial(i,j,1);
    end
end
for i=1:nx
    for j=1:nx
        T2(i,j) = temp_second(i,j,1);
    end
end
for i=1:nx
    for j=1:nx
        T3(i,j) = temp(i,j,1);
    end
end

figure()
contourf(X,Y,T1,10)
title("t=0.15s, z=0")
xlabel('x')
ylabel('y')
colorbar

figure()
contourf(X,Y,T2,10)
title("t=2.5s, z=0")
xlabel('x')
ylabel('y')
colorbar

figure()
contourf(X,Y,T3,10)
title("Steady, z=0")
xlabel('x')
ylabel('y')
colorbar

% plots at y=0
for i=1:nx
    for j=1:nx
        T1(i,j) = temp_initial(i,1,j);
    end
end
for i=1:nx
    for j=1:nx
        T2(i,j) = temp_second(i,1,j);
    end
end
for i=1:nx
    for j=1:nx
        T3(i,j) = temp(i,1,j);
    end
end

figure()
contourf(X,Y,T1,10)
title("t=0.15s,y=0")
xlabel('x')
ylabel('z')
colorbar

figure()
contourf(X,Y,T2,10)
title("t=2.5s, y=0")
xlabel('x')
ylabel('z')
colorbar

figure()
contourf(X,Y,T3,10)
title("Steady, y=0")
xlabel('x')
ylabel('z')
colorbar


% plots at x=0
for i=1:nx
    for j=1:nx
        T1(i,j) = temp_initial(1,i,j);
    end
end
for i=1:nx
    for j=1:nx
        T2(i,j) = temp_second(1,i,j);
    end
end
for i=1:nx
    for j=1:nx
        T3(i,j) = temp(1,i,j);
    end
end

figure()
contourf(X,Y,T1,10)
title("t=0.15s,x=0")
xlabel('y')
ylabel('z')
colorbar

figure()
contourf(X,Y,T2,10)
title("tau=2.5s, x*=0")
xlabel('y*')
ylabel('z*')
colorbar

figure()
contourf(X,Y,T3,10)
title("Steady, x*=0")
xlabel('y*')
ylabel('z*')
colorbar

Theta1=zeros(nx,1);
Theta2=zeros(nx,1);
Theta3=zeros(nx,1);

% Theta vs y* at centerline x* = z* = 0
for i=1:nx
    Theta1(i)=temp_initial(1,i,1);
    Theta2(i)=temp_second(1,i,1);
    Theta3(i)=temp(1,i,1);
end

figure()
plot(y,Theta1,'r',y,Theta2,'b',y,Theta3,'g')
title("T vs y at centerline x = z = 0")
xlabel("y*")
ylabel("Theta")
legend('t=0.15 s','t=2.5 s','steady')


% T vs y* at x* = z* = L/2
for i=1:nx
    Theta1(i)=temp_initial((nx+1)/2,i,(nx+1)/2);
    Theta2(i)=temp_second((nx+1)/2,i,(nx+1)/2);
    Theta3(i)=temp((nx+1)/2,i,(nx+1)/2);
end

figure()
plot(y,Theta1,'r',y,Theta2,'b',y,Theta3,'g')
title("T vs y at x = z = L/2")
xlabel("y*")
ylabel("Theta")
legend('t=0.15 s','t=2.5 s','steady')
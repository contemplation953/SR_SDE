clc;clear;

dt = 0.01;
tspan=dt:dt:1000;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1));

dof_N=10;
[t,x]=ode45(@(t,x) diff_fun(t,dof_N,x),tspan,rand(dof_N*2,1),options);
figure(1)
plot(x(:,1))


function dy=diff_fun(t,N,yn)
    dy=zeros(length(yn),1);
    q=yn(1:N);
    p=yn(N+1:end);

    for i=1:N
        fi=0.01*q(i).*q(i)+0.02;
        if i==1
            gi=q(1)+q(1)-q(2)+(q(1)-q(2)).*(q(1)-q(2)).*(q(1)-q(2));
        elseif i==N
            gi=q(N)+q(N)-q(N-1)+(q(N)-q(N-1)).*(q(N)-q(N-1)).*(q(N)-q(N-1));
        else
            gi=q(i)+q(i)-q(i-1)+(q(i)-q(i-1)).*(q(i)-q(i-1)).*(q(i)-q(i-1)) ...
                +q(i)-q(i+1)+(q(i)-q(i+1)).*(q(i)-q(i+1)).*(q(i)-q(i+1));
        end
        dy(i)=p(i);
        dy(N+i)=-fi.*p(i)-gi;
    end
end
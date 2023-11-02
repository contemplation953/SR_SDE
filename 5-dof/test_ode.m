clc;clear;

dt = 0.01;
tspan=dt:dt:1000;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1));

[t,x]=ode45(@(t,x) fun(t,x),tspan,ones(10,1),options);
figure(1)
plot(x(:,1))



function dy = fun(t,yn)
    dy=zeros(length(yn),1);
    q1 = yn(1);
    q2 = yn(2);
    q3 = yn(3);
    q4 = yn(4);
    q5 = yn(5);

    p1 = yn(6);
    p2 = yn(7);
    p3 = yn(8);
    p4 = yn(9);
    p5 = yn(10);

    
    dy(1)=p1;
    dy(2)=p2;
    dy(3)=p3;
    dy(4)=p4;
    dy(5)=p5;
    dy(6)=-(0.01.*q1.*q1+0.02).*p1-(1*1.*q1+1*(q1-q2)+1.*(q1-q2).*(q1-q2).*(q1-q2));
    dy(7)=-(0.01.*q2.*q2+0.02).*p2-(1*1.*q2+1*(q2-q1)+1.*(q2-q1).*(q2-q1).*(q2-q1)+1.*(q2-q3)+1.*(q2-q3).*(q2-q3).*(q2-q3));
    dy(8)=-(0.01.*q3.*q3+0.02).*p3-(1*1.*q3+1*(q3-q2)+1.*(q3-q2).*(q3-q2).*(q3-q2)+1.*(q3-q4)+1.*(q3-q4).*(q3-q4).*(q3-q4));
    dy(9)=-(0.01.*q4.*q4+0.02).*p4-(1*1.*q4+1*(q4-q3)+1.*(q4-q3).*(q4-q3).*(q4-q3)+1.*(q4-q5)+1.*(q4-q5).*(q4-q5).*(q4-q5));
    dy(10)=-(0.01.*q5.*q5+0.02).*p5-(1*1.*q5+1*(q5-q4)+1.*(q5-q4).*(q5-q4).*(q5-q4));
end
clc;clear;


root_path=pwd();

tic;

pm=2;
ps=5;
x_list=zeros(61,pm+ps+2);

i=1;
data=load(sprintf('%s/5-dof/data/H_%d.mat',root_path,i));

h=0.01;
H_list1=0.01:h:25-h;


moment1=data.moment1(1:length(H_list1));
moment2=sqrt(data.moment2(1:length(H_list1))); 

coeff=init(H_list1,moment1,moment2,pm,ps);
H_list=0.01:h:20-h;
y1=polyval(coeff(1:pm+1),H_list);
y2=polyval(coeff(pm+2:end),H_list);
p_data=data.p_data(1:length(H_list));


% options = optimset('MaxFunEvals',1e9,'MaxIter',1e3,'TolFun',1e-4);
options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',1e9,'MaxIter',1e3,'TolFun ',1e-2);
fun=@(coeff) ssr_cost(coeff,pm,ps,p_data,H_list,y1,y2);

[x,val]=fminsearch(fun,coeff,options);

toc;

function coeff=init(pst,m1_data,m2_data,pm,ps)
    coeff1=polyfit(pst,m1_data,pm);
    coeff2=polyfit(pst,m2_data,ps);
    coeff=[coeff1,coeff2];
end


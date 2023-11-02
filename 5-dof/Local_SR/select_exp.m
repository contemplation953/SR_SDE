clc;clear;

root_path=pwd();
res=zeros(5,5);
oned=zeros(1,25);
index=1;
ssr_path=sprintf('%s/5-dof/ssr',root_path);
for i=1:5
    for j=1:5
        data=load(sprintf('%s/SSR_%d_%d.mat',ssr_path,i,j));
        res(i,j)=data.res;
        oned(index)=data.res;
        index=index+1;
    end
end


figure(1);
bar3(res,0.5,'grouped');
ylabel('p(m)');
zlabel('Loss_{pq}');
legend('p(\sigma)=1','p(\sigma)=2','p(\sigma)=3','p(\sigma)=4','p(\sigma)=5');

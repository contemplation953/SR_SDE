function root_list = solve_polynomial(f,boundery,n,eps)
% 弦截法
%
% 输入：
% f: 方程对应函数
% x0: 迭代初值
% eps: 精度
%
% 输出：
% x: 根
% X: 根序列
% iter: 迭代次数

if nargin < 4
    eps = 1e-8;
end

a=boundery(1);
b=boundery(2);

front = a;
interval_list = [];
interval = a : 1e-2: b;
for i = 1:length(interval)-1
    if f(interval(i)) == 0
        interval_list = [interval_list;[front,interval(i+1)]];
        front = interval(i+1);
    elseif f(front)*f(interval(i)) > 0
        front = interval(i);
    else
        interval_list = [interval_list;[front,interval(i)]];
        front = interval(i);
    end
end

assert(length(interval_list) == n)
root_list = zeros(1,n);
for i = 1 : n
    a=interval_list(i,1);
    b=interval_list(i,2);
    root_list(i)=Bisection(f,a,b,eps);
end

end

function x = Bisection(f,a,b,eps)
% 二分法求解非线性方程
%
% 输入：
% f: 非线性方程对应函数
% [a,b]: 求根区间
% eps: 精度
%
% 输出：
% x: 根
% [A,B]: 求根区间序列
% M: 根序列
% iter: 迭代次数

if nargin < 4
    eps = 1e-6;
end

A = [];
B = [];
X = [];

iter = 0;

while abs(a-b) > eps
    iter = iter+1;
    x = (a+b)/2;
    X = [X;x];
    A = [A;a];
    B = [B;b];
    if f(x) == 0
        break
    else
        if f(a)*f(x)<0
            b = x;
        else
            a = x;
        end
    end
end
end
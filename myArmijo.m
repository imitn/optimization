function [mk,alpha,fk,newxk,newfk] =myArmijo(xk,dk,fun,gfun)

%功能：Armijo准则非精确线搜索
%输入：xk是起始变量，
%     dk是下降方向，
%     fun是目标函数
%     gfun是目标函数的导数
%输出：mk是满足不等式的最小非负整数，
%     alpha是步长，
%     fk是初始变量的函数值，
%     newxk是当前迭代变量
%     newfk是当前迭代变量的目标函数值

%例子：>> fun=@(x) 100*(x(1)^2-x(2))^2+(x(1)-1)^2;
%     >> gfun=@(x) [400*x(1)*(x(1)^2-x(2))+2*(x(1)-1);-200*(x(1)^2-x(2))]
%     >> xk=[-1,1]';
%     >> dk=[1,-2]';
%     >> [mk,alpha,fk,newxk,newfk] =myArmijo(xk,dk,fun,gfun);

%如果运行出现错误，matlab会自动停在出错的那行，并且保存所有相关变量
dbstop if error

beta=0.5;
sigma=0.2;
m=0;
maxm=20;

while m<maxm
    if fun(xk+beta^m*dk) <= (fun(xk)+sigma*beta^m*gfun(xk)'*dk)
        mk=m;
        break;
    end
    m=m+1;
end
alpha=beta^mk;
newxk=xk+alpha*dk;
fk=fun(xk);
newfk=fun(newxk);

end

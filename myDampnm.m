function [k,x,val] = myDampnm(x0,epsilon,fun,gfun,Hess)

%功能：阻尼牛顿法求解无约束优化问题：min f(x)
%输入：x0是初始点，
%     epsilon为容许误差，
%     fun是目标函数
%     gfun是梯度
%     Hess是Hesse阵
%输出：k是迭代次数，
%     x和val分别是近似最优点和最优值，

%例子：>> fun=@(x) 4*(x(1)^2-x(2))^2+3*(x(1)-1)^2;
%     >> gfun=@(x) [16*x(1)*(x(1)^2-x(2))+6*(x(1)-1);-8*(x(1)^2-x(2))];
%     >> Hess=@(x) [48*x(1)^2-16*x(2)+6,-16*x(1);-16*x(1),8];
%     >> x0=[-1.2,1.0]';
%     >> epsilon=1e-5';
%     >> [k,x,val] = myDampnm(x0,epsilon,fun,gfun,Hess);

%如果运行出现错误，matlab会自动停在出错的那行，并且保存所有相关变量
dbstop if error

maxk=5000;
beta=0.5;
sigma=0.4;

k=0;
while k<maxk
    
    gk=gfun(x0);%计算梯度
    Gk=Hess(x0);%计算Hesse阵
    dk=-Gk\gk;%解方程组 Gk*dk=-gk，计算搜索方向（牛顿方向）
    
    %检验终止条件，梯度变化太小就结束
    if norm(dk)<epsilon
        break
    end
    
    %使用Armijo搜索求步长
    m=0;
    mk=0;
    while m<20
       if fun(x0+beta^m*dk) <= (fun(x0)+sigma*beta^m*gk'*dk)
          mk=m;
          break;
       end
       m=m+1;
    end
    x0=x0+beta^mk*dk;
    k=k+1;    
end
x=x0;
val=fun(x0);

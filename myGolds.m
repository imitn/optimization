function [i,s,phis,ds,dphi,G] = myGolds(phi,a,b,epsilon,delta)

%功能：黄金分割法精确线性搜索
%输入：phi是目标函数，
%     a，b是搜索区间的两个端点，
%     epsilon和delta分别是自变量和函数值的容许误差
%输出：i是迭代次数，
%     s和phis分别是近似极小值点和极小值，
%     ds和dphi分别是s和phis的误差限，
%     G是i乘4大小的矩阵，其第i行分别是a，p，q，b的第i次迭代值[ai,pi,qi,bi]

%例子：>> phi=@(x) 3*x^2-2*tan(x)
%     >>function [i,s,phis,ds,dphi,G] = myGolds(phi,0,1,1e-4,1e-5)

%如果运行出现错误，matlab会自动停在出错的那行，并且保存所有相关变量
dbstop if error

t=(sqrt(5)-1)/2;
h=b-a;

phia=feval(phi,a);
phib=feval(phi,b);

p=a+(1-t)*h;
q=a+t*h;
phip=feval(phi,p);
phiq=feval(phi,q);

i=1;
G(i,:)=[a,p,q,b];

while (abs(phib-phia)>delta) || (h>epsilon)
    if phip<=phib %计算左试探点
        %右侧区间向左收缩，开始
        b=q;
        phib=phiq;
        q=p;
        phiq=phip;
        %结束
        %计算试探点
        h=b-a;
        p=a+(1-t)*h;
        phip=feval(phi,p);     
    else %计算右试探点
        %右侧区间向右收缩，开始
        a=p;
        phia=phip;
        p=q;
        phip=phiq;
        %结束
        %计算试探点
        h=b-a;
        q=a+t*h;
        phiq=feval(phi,q);           
    end
    i=i+1;
    G(i,:)=[a,p,q,b];
end

if phip<=phiq
    s=p;
    phis=phip;
else
    s=q;
    phis=phiq;
end

ds=abs(b-a);
dphi=abs(phib-phia);

end

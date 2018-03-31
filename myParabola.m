function [i,s,phis,ds,dphi,S] = myParabola( phi,a,b,epsilon,delta )

%功能：抛物线法精确线性搜索
%输入：phi是目标函数，
%     a，b是搜索区间的两个端点，
%     epsilon和delta分别是自变量和函数值的容许误差
%输出：i是迭代次数，
%     s和phis分别是近似极小值点和极小值，
%     ds和dphi分别是|s-s1|和|phi(s)-phi(s1)|，
%     S是向量，其第i个分量是i次迭代值s0

%例子：>> phi=@(x) 3*x^2-2*tan(x)
%     >>function [i,s,phis,ds,dphi,S] = myParabola(phi,0,1,1e-4,1e-5)

%如果运行出现错误，matlab会自动停在出错的那行，并且保存所有相关变量
dbstop if error

s0=a;
i=1;
S(i)=s0;
maxi=30;
maxj=20;
big=1e6;
err=1;
cond=0;
h=1;
ds=0.00001;
if abs(s0)>1e4
    h=abs(s0)*(1e-4);
end

while i<maxi && err>delta && cond~=5
    
    f1=(feval(phi,s0+ds)-feval(phi,s0-ds))/(2*ds);
    
    if f1>0
        h=abs(h);
    end
    
    s1=s0+h;
    s2=s0+2*h;
    bars=s0;
    
    phi0=feval(phi,s0);
    phi1=feval(phi,s1);
    phi2=feval(phi,s2);
    barphi=phi0;
    cond=0;
    
    %确定h使得phi1<phi0且phi1<phi2
    j=0;
    while j<maxj && abs(h)>epsilon && cond==0
        if phi0<=phi1
            s2=s1;
            phi2=phi1;
            h=0.5*h;
            s1=s0+h;
            phi1=feval(phi,s1);
        elseif phi2<phi1
             s1=s2;
             phi1=phi2;
             h=2*h;
             s2=s0+2*h;
             phi2=feval(phi,s2);
        else
           cond=-1; 
        end   

        j=j+1;
        if abs(h)>big || abs(s0)>big
           cond=5; 
        end
    end
    
    if cond==5
        bars=s1;
        barphi=feval(phi,s1);
    else
        %二次插值求phis
        d=2*(2*phi1-phi0-phi2);
        if d<0
            barh=h*(4*phi1-3*phi0-phi2)/d;
        else
            barh=h/3;
            cond=4;
        end
        bars=s0+barh;
        barphi=feval(phi,bars);
        h=abs(h);
        h0=abs(barh);
        h1=abs(barh-h);
        h2=abs(barh-2*h);
        
        %确定下一次迭代的h值
        if h0<h
            h=h0;
        end
        if h1<h
            h=h1;
        end
        if h2<h
            h=h2;
        end   
        if h==0
            h=barh;
        end
        if h<epsilon
            cond=1;
        end
        if abs(h)>big || abs(bars)>big
            cond=5;
        end
        err=abs(barphi-phi1);
        s0=bars;
        i=i+1;
        S(i)=s0;
    end
    if cond==2 && h<epsilon
        cond=3;
    end
end

s=s0;
phis=feval(phi,s);
ds=abs(s-s1);
dphi=err;

end

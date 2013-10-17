close all; clc;

rng(0);       % regularize randomization for testing

n=41; %原子层数
m=41; %长10nm
h=0.25;  %单位nm
time=5*10^-8;
deltat=10^-9; %时间步长s
Ei=zeros(1,m);   %单位v/nm
Ec=zeros(1,m);   %单位v/nm
Ea=1;   %势垒高度，单位eV
kb=1.38*10^-23; %boltzmann常数
T=300;  %温度，单位k
qe=1.6*10^-19;
t0=10^-13;  %氧离子本征频率
pg=zeros(1,m);  %每列产生率
voltage=2;    %V
ratio=100;  %氧空位与绝缘层的电阻率之比
velocity=zeros(1,m);    %每列氧离子迁移速率
vo= false(n,m);  %氧空位矩阵，1代表该位置有空位
io= false(n,m);  %氧离子矩阵，1代表该位置有离子
ganma=3;  %e*nm
miu=2.5*10^9;   %nm^2/(V*s)
history_vo={};
iocolumn=zeros(1,m);    %一列中的氧离子
for v=1.1:0.1:voltage
for t=1:1:time/deltat
    vocolumn=zeros(1,m);    %一列中的氧空位
    for i=1:n
        for j=1:m
            if(vo(i,j)==1)
                vocolumn(1,j)=vocolumn(1,j)+1;
            end
        end
    end
    for i=1:m
       Ei(1,i)=v/((n-1)*h-vocolumn(1,i)*h)/(1+vocolumn(1,i)*h/ratio/((n-1)*h-vocolumn(1,i)*h));
       pg(1,i)=deltat/t0*exp(-(Ea-ganma*Ei(1,i))/(kb*T/qe));
       velocity(1,i)=miu*Ei(1,i);
    end
    for i=1:n
        for j=1:m
            if(io(i,j)==1)
                length=floor(velocity(1,j)*deltat/h);
                recombination=0;
                if(i==1)
                    io(i,j)=0;
                    iocolumn(1,j)=iocolumn(1,j)+1;
                else
                    for k=i-1:-1:max(1,i-length)
                        if(vo(k,j)==1)
                            vo(k,j)=0;
                            io(i,j)=0;
                            recombination=1;
                            break;
                        end
                    end
                    if(recombination==0)
                        if(i-length>0)
                            io(i,j)=0;
                            io(i-length,j)=1;
                        else
                            io(i,j)=0;
                            iocolumn(1,j)=iocolumn(1,j)+1;
                        end
                    end
                end
            end
            r=unifrnd(0,1);
            if(vo(i,j)==0&&r<pg(1,j))
                vo(i,j)=1;
                io(i,j)=1;
            end
        end
    end
end
history_vo=[history_vo,{vo}];
end
save('distribution','history_vo');

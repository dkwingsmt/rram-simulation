clear all;close all; clc;


load_consts

time=2*10^-7;
deltat=2*10^-9; %时间步长s

Ei=zeros(1,m);   %单位V/m
Ec=zeros(1,m);   %单位V/m
pg=zeros(1,m);  %每列产生率
velocity=zeros(1,m);    %每列氧离子迁移速率
history_vo={};
iocolumn=zeros(1,m);    %一列中的氧离子

voltage=2.6;    %V
for v=1.1:0.1:voltage
vo= zeros(n,m);  %氧空位矩阵，1代表该位置有空位
io= zeros(n,m);  %氧离子矩阵，1代表该位置有离子
iocolumn=zeros(1,m);    %一列中的氧离子
for t=1:1:time/deltat

    %% Calculate Vo configuration
    vocolumn = sum(vo, 1);  % Amount of Vo for every column

    Ei = v ./ ((n-1)*h-vocolumn*h) ./ (1+vocolumn*h/ratio/((n-1)*h-vocolumn*h));
    pg = deltat/t0*exp(-(Ea-gamma*Ei)/(kb*T));
    velocity = miu*Ei;

    r=unifrnd(0,1,1,m);   
    for i=1:n
        for j=1:m
            if(io(i,j)==1)
                max_distance=floor(velocity(1,j)*deltat/h);
                recombination=0;
                if(i==1)
                    iocolumn(1,j)=iocolumn(1,j)+1;
                else
                    for k=i-1:-1:max(1,i-max_distance)
                        if(vo(k,j)==1)
                            vo(k,j)=0;
                            recombination=1;
                            break;
                        end
                    end
                    if(recombination==0)
                        if(i-max_distance>0)
                            io(i-max_distance,j)=1;
                        else
                            iocolumn(1,j)=iocolumn(1,j)+1;
                        end
                    end
                end
                io(i,j)=0;
            end
            if(vo(i,j)==0&&r(1,j)<pg(1,j))
                vo(i,j)=1;
                io(i,j)=1;
            end
        end
    end
    if(max(iocolumn)==n)
        break;
    end
end
for i=1:n
    for j=1:m
        if(vo(i,j)==1&&io(i,j)==1)
            vo(i,j)=0;
            io(i,j)=0;
        end
    end
end
history_vo=[history_vo,{vo}];
end
save('distribution','history_vo');

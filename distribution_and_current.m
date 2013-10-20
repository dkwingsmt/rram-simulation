clear all;close all; clc;


load_consts

time=2*10^-7;
deltat=2*10^-9; %ʱ�䲽��s

Ei=zeros(1,m);   %��λV/m
Ec=zeros(1,m);   %��λV/m
pg=zeros(1,m);  %ÿ�в�����
velocity=zeros(1,m);    %ÿ��������Ǩ������
history_vo={};
iocolumn=zeros(1,m);    %һ���е�������

voltage=2.6;    %V
for v=1.1:0.1:voltage
vo= zeros(n,m);  %����λ����1�����λ���п�λ
io= zeros(n,m);  %�����Ӿ���1�����λ��������
iocolumn=zeros(1,m);    %һ���е�������
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

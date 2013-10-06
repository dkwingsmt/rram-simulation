close all; clc;

rng(0);       % regularize randomization for testing

n=41; %ԭ�Ӳ���
m=41; %��10nm
h=0.25;  %��λnm
time=5*10^-8;
deltat=10^-9; %ʱ�䲽��s
Ei=zeros(1,m);   %��λv/nm
Ec=zeros(1,m);   %��λv/nm
Ea=1;   %���ݸ߶ȣ���λeV
kb=1.38*10^-23; %boltzmann����
T=300;  %�¶ȣ���λk
qe=1.6*10^-19;
t0=10^-13;  %�����ӱ���Ƶ��
pg=zeros(1,m);  %ÿ�в�����
voltage=2;    %V
ratio=100;  %����λ���Ե��ĵ�����֮��
velocity=zeros(1,m);    %ÿ��������Ǩ������
vo= false(n,m);  %����λ����1�����λ���п�λ
io= false(n,m);  %�����Ӿ���1�����λ��������
ganma=3;  %e*nm
miu=2.5*10^9;   %nm^2/(V*s)
history_vo={};
iocolumn=zeros(1,m);    %һ���е�������
for v=1.1:0.1:voltage
for t=1:1:time/deltat

    %% Calculate Vo configuration
    vocolumn = sum(vo, 1);  % Amount of Vo for every column

    Ei = v ./ ((n-1)*h-vocolumn*h) ./ (1+vocolumn*h/ratio/((n-1)*h-vocolumn*h));
    pg = deltat/t0*exp(-(Ea-ganma*Ei)/(kb*T/qe));
    velocity = miu*Ei;

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
        end
    end

    r = unifrnd(0, 1, n, m);
    do_generate = (vo(i,j)==0) & (r < repmat(pg, [n, 1])); 
    vo(do_generate) = 1;
    io(do_generate) = 1;

end
history_vo=[history_vo,{vo}];
end
save('distribution','history_vo');

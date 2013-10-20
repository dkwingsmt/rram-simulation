close all; clear all; clc;

load distribution
load_consts

plotnum = length(history_vo);
for plot_id = 1:plotnum
        figure(1);
        vo=history_vo{plot_id};
        x=1:m+1;
        y=1:n+1;
        M=meshgrid(x,y);
        N=meshgrid(y,x);
        plot(x,N,'k');
        hold on;
        plot(M,y,'k');
        for i=1:n
            for j=1:m
                if(vo(i,j)==1)
                    x=[j,j+1,j+1,j];
                    y=[i,i,i+1,i+1];
                    fill(x,y,'r');
                end
                hold on;
            end
        end
        hold off;
end
function [D,IXY,IP,GX]=G2D(G)
% 生成任意两个近邻栅格间的距离,相邻的无障碍的有距离,相邻的有障碍的为0

G=flip(G);
% 可行的区域进行编号
[km,kn]=size(G);
GX=1-G;
L=0;
for m=1:km
    for n=1:kn
        if abs(GX(m,n)-1)<1e-8
            L=L+1;
            GX(m,n)=L;
        end
    end
end
DN=L;
% 周围节点间距(左上,上,右上,左,右,左下,下,右下)
d=[sqrt(2),1,sqrt(2),1,1,sqrt(2),1,sqrt(2)];
D=zeros(DN,length(d)); 
% 节点坐标
for i=1:DN
    % 该点的位置
    [Iy,Ix]=find(abs(i-GX)<1e-8);
    % 该点的坐标
    IXY(i,:)=[Ix,Iy];
    % 近邻栅格的编号(左上,上,右上,左,右,左下,下,右下)
    IS=[Ix-1,Iy-1;Ix,Iy-1;Ix+1,Iy-1;Ix-1,Iy;Ix+1,Iy;...
        Ix-1,Iy+1;Ix,Iy+1;Ix+1,Iy+1];
    % 求距离和节点
    for sn=1:size(IS,1)
        if (IS(sn,1)<=0.5 | IS(sn,2)<=0.5 | IS(sn,1)>kn+0.5 | IS(sn,2)>km+0.5)
            D(i,sn)=0;
            IP(i,sn)=0;
        else
            IP(i,sn)=GX(IS(sn,2),IS(sn,1));
            D(i,sn)=d(sn)*(IP(i,sn)>0.5);
        end
    end
end
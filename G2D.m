function [D,IXY,IP,GX]=G2D(G)
% ����������������դ���ľ���,���ڵ����ϰ����о���,���ڵ����ϰ���Ϊ0

G=flip(G);
% ���е�������б��
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
% ��Χ�ڵ���(����,��,����,��,��,����,��,����)
d=[sqrt(2),1,sqrt(2),1,1,sqrt(2),1,sqrt(2)];
D=zeros(DN,length(d)); 
% �ڵ�����
for i=1:DN
    % �õ��λ��
    [Iy,Ix]=find(abs(i-GX)<1e-8);
    % �õ������
    IXY(i,:)=[Ix,Iy];
    % ����դ��ı��(����,��,����,��,��,����,��,����)
    IS=[Ix-1,Iy-1;Ix,Iy-1;Ix+1,Iy-1;Ix-1,Iy;Ix+1,Iy;...
        Ix-1,Iy+1;Ix,Iy+1;Ix+1,Iy+1];
    % �����ͽڵ�
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
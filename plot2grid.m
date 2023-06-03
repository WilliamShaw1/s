function plot2grid(G,GN,ROUT)
% 画出栅格和比较优的爬行路径

if nargin<=2;   ROUT=[];   end
if nargin<=1;   GN=[];      end
[Gy,Gx]=size(G);

% 画出栅格
b=flip(G);
b(end+1,end+1)=0;
% 创建颜色
colormap([1,1,1;0.2,0.2,0.2]);  
% 赋予栅格颜色
pcolor(b); 
hold on;
if Gy*Gx>(10000+0.5)
    shading interp;
end
% 设置坐标
XL=round(linspace(1,Gx,10));
YL=round(linspace(1,Gy,10));
set(gca,'XTick',XL+0.5,'XTickLabel',XL,'YTick',YL+0.5,'YTickLabel',YL);
xlabel('坐标 x');
ylabel('坐标 y');

% 画出最短爬行路径
if (length(GN)>0.5 & length(ROUT)>0.5)
    for m=1:length(ROUT)
        [Iy,Ix]=find(abs(ROUT(m)-GN)<1e-8);
        RXY(m,:)=[Ix,Iy]+0.5;
    end
    plot(RXY(:,1),RXY(:,2),'b');
    hold on;
    % 画出起点和终点
    plot(RXY(1,1),RXY(1,2),'r*');
    hold on;
    plot(RXY(end,1),RXY(end,2),'r*');
end
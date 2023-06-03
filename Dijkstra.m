function [D,S,R]=Dijkstra(W)
% 改进的Dijkstra算法,求解所有两点之间的最短距离

% 求出节点数
N=size(W,1);
D=zeros(N);
% 路径初始化
R=cell(N,N);
for n=1:N
    % 初始时,每一点到各点的距离
    DA(1:N)=max(max(W));
    % 该点到自己的距离为0
    DA(n)=0;
    % 已找到最短路径的节点集合
    SA=zeros(1,N);
    SA(1)=n;
    % 未找到最短路径的节点集合
    T=setdiff([1:N],SA);
    % 最短距离的路径开辟存储空间
    RA=cell(1,N);
    RA{n}=[n,n];
    % 顶点之前的一个顶点
    index=zeros(1,N);
    index(n)=n;
    % 参数初始化
    temp=n;
    k=1;
    ps=SA(SA>=1);
    % 更新l(v),同时记录顶点顺序和顶点索引
    while length(ps)<N
        k=k+1;
        DD=[DA(T);DA(temp)+W(temp,T)];
        index(T)=(DD(1,:)<=DD(2,:)).*index(T)+(DD(1,:)>DD(2,:)).*temp;
        DA(T)=min(DD);
        [~,pd]=min(DA(T));
        temp=T(pd(1));
        SA(k)=temp;
        T=setdiff(T,temp);
        ps=SA(SA>=1);
    end
    % 源点到各点最短距离的路径
    for m=setdiff([1:N],n)
        IN=index(m);
        if abs(IN)<1e-8
            RA{m}=NaN;
        else
            route=[IN,m];
            while abs(IN-n)>1e-8
                IN=index(IN);
                route=[IN,route];
            end
            RA{m}=route;
        end
    end
    D(n,:)=DA;
    S(n,:)=SA;
    R(n,:)=RA;
end
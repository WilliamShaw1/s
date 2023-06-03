function [Rbest,SRN,SR,SL,Lbest,Lave]=ant1Dad(CT,W,M,IT,alpha,beta,rho,Q)
% 一维蚁群算法,用于旅行商问题(TSP)优化,从起点出发再回到起点

% 城市名
CT=reshape(CT,length(CT),1);
% 城市数目
N=size(CT,1);
% 生成任意两城市之间的距离
D=W;
% 对角线上的元素
for n=1:N
    D(n,n)=1e-8;
end

% 初始化设置
Eta=1./D;                               %启发函数
Tau=ones(N);                            %信息素矩阵
Table=zeros(M,N);                       %路径记录表
Rbest=zeros(IT,N);                      %各代最佳路径
Lbest=zeros(IT,1);                      %各代最佳路径的长度
Lave=zeros(IT,1);                       %各代路径的平均长度
iter=1;                                 %迭代次数初值

% 迭代寻找最佳路径
while iter<=IT
    % 随机产生各个蚂蚁的起点城市
    st=randi([1,N],M,1);
    Table(:,1)=st;
    % 构建解空间
    CIN=[1:N];
    % 逐个蚂蚁路径选择
    for m=1:M
        % 逐个城市路径选择
        for n=2:N
            %已访问的城市集合(禁忌表)
            tabu=Table(m,1:(n-1));          
            %待访问的城市集合
            allow=setdiff(CIN,tabu); 
            P=allow;
            % 计算城市间转移概率
            for k=1:length(allow)
                P(k)=Tau(tabu(end),allow(k))^alpha*Eta(tabu(end),allow(k))^beta;
            end
            P=P/sum(P);
            % 轮盘赌法选择下一个访问城市
            Pc=cumsum(P);
            TAI=find(Pc>=rand);
            TAR=allow(TAI(1));
            Table(m,n)=TAR;
        end
    end
    % 计算各个蚂蚁的路径距离
    Len=zeros(M,1);
    for m=1:M
        Route=[Table(m,:),Table(m,1)];
        ch1=Route(1:end-1);
        ch2=Route(2:end);
        d=D(ch1,ch2);
        Len(m)=sum(diag(d));
    end
    % 计算最短路径距离及平均距离
    if abs(iter-1)<1e-8
        [minL,minIN]=min(Len);
        Lbest(iter)=minL;
        Lave(iter)=mean(Len);
        Rbest(iter,:)=Table(minIN,:);
    else
        [minL,minIN]=min(Len);
        Lbest(iter)=min(Lbest(iter-1),minL);
        Lave(iter)=mean(Len);
        if abs(Lbest(iter)-minL)<1e-8
            Rbest(iter,:)=Table(minIN,:);
        else
            Rbest(iter,:)=Rbest((iter-1),:);
        end
    end
    % 更新信息素
    Delta_Tau=zeros(N);
    % 逐个蚂蚁计算
    for m=1:M
        % 逐个城市计算
        for n=1:(N-1)
            i=Table(m,n);
            j=Table(m,n+1);
            Delta_Tau(i,j)=Delta_Tau(i,j)+Q/Len(m);
        end
        i=Table(m,N);
        j=Table(m,1);
        Delta_Tau(i,j)=Delta_Tau(i,j)+Q/Len(m);
    end
    Tau=(1-rho)*Tau+Delta_Tau;
    % 迭代次数加1,清空路径记录表
    iter=iter+1;
    Table=zeros(M,N);
end

% 结果显示
[SL,index]=min(Lbest);      %最短距离
SR=[Rbest(index,:),Rbest(index,1)];          %最短路径
SRN=CT(SR,1);
SRN=reshape(SRN,1,length(SRN));

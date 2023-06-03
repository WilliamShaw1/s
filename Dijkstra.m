function [D,S,R]=Dijkstra(W)
% �Ľ���Dijkstra�㷨,�����������֮�����̾���

% ����ڵ���
N=size(W,1);
D=zeros(N);
% ·����ʼ��
R=cell(N,N);
for n=1:N
    % ��ʼʱ,ÿһ�㵽����ľ���
    DA(1:N)=max(max(W));
    % �õ㵽�Լ��ľ���Ϊ0
    DA(n)=0;
    % ���ҵ����·���Ľڵ㼯��
    SA=zeros(1,N);
    SA(1)=n;
    % δ�ҵ����·���Ľڵ㼯��
    T=setdiff([1:N],SA);
    % ��̾����·�����ٴ洢�ռ�
    RA=cell(1,N);
    RA{n}=[n,n];
    % ����֮ǰ��һ������
    index=zeros(1,N);
    index(n)=n;
    % ������ʼ��
    temp=n;
    k=1;
    ps=SA(SA>=1);
    % ����l(v),ͬʱ��¼����˳��Ͷ�������
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
    % Դ�㵽������̾����·��
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
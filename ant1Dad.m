function [Rbest,SRN,SR,SL,Lbest,Lave]=ant1Dad(CT,W,M,IT,alpha,beta,rho,Q)
% һά��Ⱥ�㷨,��������������(TSP)�Ż�,���������ٻص����

% ������
CT=reshape(CT,length(CT),1);
% ������Ŀ
N=size(CT,1);
% ��������������֮��ľ���
D=W;
% �Խ����ϵ�Ԫ��
for n=1:N
    D(n,n)=1e-8;
end

% ��ʼ������
Eta=1./D;                               %��������
Tau=ones(N);                            %��Ϣ�ؾ���
Table=zeros(M,N);                       %·����¼��
Rbest=zeros(IT,N);                      %�������·��
Lbest=zeros(IT,1);                      %�������·���ĳ���
Lave=zeros(IT,1);                       %����·����ƽ������
iter=1;                                 %����������ֵ

% ����Ѱ�����·��
while iter<=IT
    % ��������������ϵ�������
    st=randi([1,N],M,1);
    Table(:,1)=st;
    % ������ռ�
    CIN=[1:N];
    % �������·��ѡ��
    for m=1:M
        % �������·��ѡ��
        for n=2:N
            %�ѷ��ʵĳ��м���(���ɱ�)
            tabu=Table(m,1:(n-1));          
            %�����ʵĳ��м���
            allow=setdiff(CIN,tabu); 
            P=allow;
            % ������м�ת�Ƹ���
            for k=1:length(allow)
                P(k)=Tau(tabu(end),allow(k))^alpha*Eta(tabu(end),allow(k))^beta;
            end
            P=P/sum(P);
            % ���̶ķ�ѡ����һ�����ʳ���
            Pc=cumsum(P);
            TAI=find(Pc>=rand);
            TAR=allow(TAI(1));
            Table(m,n)=TAR;
        end
    end
    % ����������ϵ�·������
    Len=zeros(M,1);
    for m=1:M
        Route=[Table(m,:),Table(m,1)];
        ch1=Route(1:end-1);
        ch2=Route(2:end);
        d=D(ch1,ch2);
        Len(m)=sum(diag(d));
    end
    % �������·�����뼰ƽ������
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
    % ������Ϣ��
    Delta_Tau=zeros(N);
    % ������ϼ���
    for m=1:M
        % ������м���
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
    % ����������1,���·����¼��
    iter=iter+1;
    Table=zeros(M,N);
end

% �����ʾ
[SL,index]=min(Lbest);      %��̾���
SR=[Rbest(index,:),Rbest(index,1)];          %���·��
SRN=CT(SR,1);
SRN=reshape(SRN,1,length(SRN));

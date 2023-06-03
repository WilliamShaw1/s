function [RSN,RS,DS,GCLE,Lav]=GATSPad(Chrom,CT,W,IT,Ps,Pc,Pm)
% �Ŵ��㷨���TSP����,���������ٻص����

CT=reshape(CT,length(CT),1);
N=size(CT,1);       % ������Ŀ
% Ϊÿһ�����ֵ���ٿռ�(GCL��1:N+1��ʾ·�����,GCL��N+2��ʾ��·���ĳ���)
GCL=zeros(IT,N+2);
% ��Ⱥ��Ŀ
M=size(Chrom,1);
% �ݻ�����
gen=1;
while gen<=IT
    % ������Ӧ�ȣ�·�߳��ȣ�
    ObjV=PathLength(Chrom,W);
    Lav(gen,1)=mean(ObjV);
    % �洢ÿһ������Сֵ
    [~,PO]=min(ObjV);
    GCL(gen,:)=[Chrom(PO,:),Chrom(PO,1),ObjV(PO)];
    % ѡ�����
    CH1=Select(Chrom,ObjV,Ps);
    % �������
    CH2=GAcross(Chrom,Pc);
    % �������
    CH3=Mutate(Chrom,Pm);
    % ��ת����
    CH4=Reverse(Chrom,W);
    % ��λ����
    CH5=Gression(Chrom,W);
    % �������
    CH6=Recombin(Chrom,W);
    % �ز����Ӵ�������Ⱥ
    CH=[CH1;CH2;CH3;CH4;CH5;CH6];
    % ɾ���ظ��Ļ������
    MM=[];
    for kk=2:size(CH,1)
        ds=abs(CH(1:kk-1,:)-repmat(CH(kk,:),kk-1,1));
        DS=sum(ds,2);
        po=find(DS<1e-8);
        if length(po)>0.5
            MM=[MM,kk];
        end
    end
    CH(MM,:)=[];
    % ������Ӧ�ȣ�·�߳��ȣ�
    ObjV=PathLength(CH,W);
    % ����ǰM���������
    [ObjV,IN]=sort(ObjV);
    if length(ObjV)>=M
        rr=IN([1:M]);
        Chrom=CH(rr,:);
    else
        Chrom=CH;
    end
    % ���µ�������
    gen=gen+1 ;
end
% �ҵ����Ž�
ObjV=PathLength(Chrom,W);  %����·�߳���
[minObjV,minInd]=min(ObjV);
% ���Ž��
RS=[Chrom(minInd(1),:),Chrom(minInd(1),1)];
RSN=CT(RS,1)';
DS=minObjV;
% ��������ͼ
GCLE=GCL(:,end);

function CH=Select(Chrom,ObjV,Gg)
% ѡ�����,һ��ѡ��ǰ20%�����������б���
% Ĭ�ϵ�������
if nargin<=2.5 | isempty(Gg), Gg=0.2; end
% ��Ⱥ��Ŀ
N=size(Chrom,1);
% ���������Ŀ
NS=max(floor(N*Gg+0.5),2);
% ������
[~,OI]=sortrows(ObjV,1);
po=OI(1:NS);
% ��ѡ��ĸ���
CH=Chrom(po,:);
return

function CH=Reverse(Chrom,D)
% ��ת����,���и����������ת����,
% ��Ⱥ��ĿN�ͻ��򳤶�L
[N,L]=size(Chrom);
% ����·������
ObjV=PathLength(Chrom,D);
% ��ת����
CH=Chrom;
for n=1:N
    % ��Ҫ��ת��һ���ַ�����ֹ����λ��
    rn=randperm(L,2);
    r1=min(rn);
    r2=max(rn);
    CH(n,r1:r2)=CH(n,r2:-1:r1);
end
% ������ת֮���·������
OV=PathLength(CH,D); 
IN=(OV>=ObjV);
% ������ת֮��Ľ��
CH(IN,:)=Chrom(IN,:);
return

function CH=Recombin(Chrom,D)
% ������������и���������������,
% ��Ⱥ��ĿN�ͻ��򳤶�L
[N,L]=size(Chrom);
% ����·������
ObjV=PathLength(Chrom,D);
% �������
CH=Chrom;
for n=1:N
    % ��Ҫ����ĵ�ǰ����
    ch=Chrom(n,:);
    % ��Ҫ�����һ���ַ�����ֹ����λ��
    rn=randperm(L,2);
    r1=min(rn);
    r2=max(rn);
    % �������
    chtemp=ch([r1:r2]);
    temp=chtemp(1);
    tp=temp;
    chre=chtemp([2:end]);
    chlen=length(chre);
    % ���ճ���֮��ľ���Զ�����в���
    while chlen>=1.5
        LEN=D(temp,chre);
        [LENmin,tempID]=min(LEN);
        tp=[tp,chre(tempID)];
        chre(tempID)=[];
        chlen=length(chre);
    end
    tp=[tp,chre];
    CH(n,r1:r2)=tp;
end
% ��������֮���·������
OV=PathLength(CH,D);
IN=(OV>=ObjV);
% ��������֮��Ľ��
CH(IN,:)=Chrom(IN,:);
return

function Len=PathLength(Chrom,D)
% ���㻷·����
% �����е�˳�򲹳ɻ�·
Chrom=[Chrom,Chrom(:,1)];
% ���е���Ŀ
N=size(D,1);
% ��Ⱥ����Ŀ
M=size(Chrom,1);
% Ϊÿ������ĳ��ȿ��ٿռ�
Len=zeros(M,1);
% ���㳤��
for m=1:M
    ch=Chrom(m,:);
    ch1=ch(1:end-1);
    ch2=ch(2:end);
    d=D(ch1,ch2);
    Len(m)=sum(diag(d));
end
return

function CH=Mutate(Chrom,Pm)
% �����λ������һ��ѡ��20%�ĸ�����в���
% Ĭ�ϵı������
if nargin<=1.5 | isempty(Pm), Pm=0.2; end
% ��Ⱥ��ĿN�ͻ��򳤶�L
[N,L]=size(Chrom);
% Ϊ����֮����¸��忪�ٴ洢�ռ�
mn=round(Pm*N);         % ��Ҫ����ĸ�����Ŀ
mC=randperm(N,mn);      % ����ĸ���
CH=zeros(mn,L);         % �������¸���
% ���б���
for n=1:mn
    % ��Ҫ����ĵ�ǰ����
    mc=mC(n);
    ch=Chrom(mc,:);
    % ��Ҫ��λ�������ַ���λ��
    rn=randperm(L,2);
    r1=rn(1);
    r2=rn(2);
    % ���廻λ�������ַ�
    Sr1=ch(r1);
    Sr2=ch(r2);
    % ���л�λ
    ch(r1)=Sr2;
    ch(r2)=Sr1;
    % ������
    CH(n,:)=ch;
end
return

function CH=Gression(Chrom,D)
% ��λ���������и����������λ����,
% ��Ⱥ��ĿN�ͻ��򳤶�L
[N,L]=size(Chrom);
% ����·������
ObjV=PathLength(Chrom,D);
% ������λ
CH=Chrom;
for n=1:N
    % ��Ҫ��λ�ĵ�ǰ����
    ch=Chrom(n,:);
    % ��Ҫ��λ��һ���ַ�����ֹ����λ��
    rn=randperm(L,2);
    r1=min(rn);
    r2=max(rn);
    % �ַ�����λ�Ĳ���
    chtemp=ch([(r1+1):r2]);
    chtemp=[chtemp,ch(r1)];
    CH(n,r1:r2)=chtemp;
end
% ��������֮���·������
OV=PathLength(CH,D);
IN=(OV>=ObjV);
% ��������֮��Ľ��
CH(IN,:)=Chrom(IN,:);
return

function CH=GAcross(Chrom,Pc)
% ���������������һ��ѡ��80%�ĸ�������в���
% Ĭ�ϵĽ������
if nargin<=1.5 | isempty(Pc), Pc=0.8; end
% ��Ⱥ��Ŀ
[N,L]=size(Chrom);
% �໥�������������
[X,Y]=meshgrid([1:N]);
XY=[X(:),Y(:)];
% ɾ������������ͬ�����
po=find(abs(XY(:,1)-XY(:,2))<1e-8);
XY(po,:)=[];
% ���Խ�����ܸ�����Ŀ
XYN=size(XY,1);
% Ϊ����֮����¸��忪�ٴ洢�ռ�
mn=round(Pc*N);         % ��Ҫ����ĸ�����Ŀ
mC=randperm(XYN,mn);    % ����ĸ�����
CH=zeros(2*mn,L);       % ����֮����¸���
% ���н���
for n=1:mn
    % ��Ҫ����ĵ�ǰ������
    mc=mC(n);
    r1=XY(mc,1);
    r2=XY(mc,2);
    ch1=Chrom(r1,:);
    ch2=Chrom(r2,:);
    % ����֮��Ľ��
    [A,B]=intercross(ch1,ch2);
    CH(2*n-1,:)=A;
    CH(2*n,:)=B;
end
return

function [an,bn]=intercross(a,b)
% ����������ӳ���
% ���е���Ŀ
L=length(a);
% ���������λ��
rn=randperm(L,2);
% ����Ĳ���(��r1��ʼ,��r2��ֹ)
if abs(rn(1)-rn(2))>0.5
    r1=min(rn);
    r2=max(rn);
else
    r1=rn(1);
    r2=L;
end
r12=r2-r1+1;
% ����֮ǰ�ı���
a0=a;
b0=b;
% �໥��ͻ�Ĳ���
ab0=intersect(a0(r1:r2),b0(r1:r2));
abua=setdiff(b0(r1:r2),ab0);
abub=setdiff(a0(r1:r2),ab0);
% ����֮���
an=a;
bn=b;
an(r1:r2)=zeros(r12,1);
bn(r1:r2)=zeros(r12,1);
% Ѱ�ҳ�ͻ�Ĳ���
for m=1:length(abua)
    po=find(abs(an-abua(m))<1e-8);
    an(po)=abub(m);
end
for m=1:length(abub)
    po=find(abs(bn-abub(m))<1e-8);
    bn(po)=abua(m);
end
% ���ɽ���֮���
an(r1:r2)=b0(r1:r2);
bn(r1:r2)=a0(r1:r2);
return
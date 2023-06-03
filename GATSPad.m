function [RSN,RS,DS,GCLE,Lav]=GATSPad(Chrom,CT,W,IT,Ps,Pc,Pm)
% 遗传算法求解TSP问题,从起点出发再回到起点

CT=reshape(CT,length(CT),1);
N=size(CT,1);       % 城市数目
% 为每一代最佳值开辟空间(GCL的1:N+1表示路径编号,GCL的N+2表示该路径的长度)
GCL=zeros(IT,N+2);
% 种群数目
M=size(Chrom,1);
% 演化迭代
gen=1;
while gen<=IT
    % 计算适应度（路线长度）
    ObjV=PathLength(Chrom,W);
    Lav(gen,1)=mean(ObjV);
    % 存储每一代的最小值
    [~,PO]=min(ObjV);
    GCL(gen,:)=[Chrom(PO,:),Chrom(PO,1),ObjV(PO)];
    % 选择操作
    CH1=Select(Chrom,ObjV,Ps);
    % 交叉操作
    CH2=GAcross(Chrom,Pc);
    % 变异操作
    CH3=Mutate(Chrom,Pm);
    % 逆转操作
    CH4=Reverse(Chrom,W);
    % 移位操作
    CH5=Gression(Chrom,W);
    % 重组操作
    CH6=Recombin(Chrom,W);
    % 重插入子代的新种群
    CH=[CH1;CH2;CH3;CH4;CH5;CH6];
    % 删除重复的基因编码
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
    % 计算适应度（路线长度）
    ObjV=PathLength(CH,W);
    % 保留前M个优秀个体
    [ObjV,IN]=sort(ObjV);
    if length(ObjV)>=M
        rr=IN([1:M]);
        Chrom=CH(rr,:);
    else
        Chrom=CH;
    end
    % 更新迭代次数
    gen=gen+1 ;
end
% 找到最优解
ObjV=PathLength(Chrom,W);  %计算路线长度
[minObjV,minInd]=min(ObjV);
% 最优结果
RS=[Chrom(minInd(1),:),Chrom(minInd(1),1)];
RSN=CT(RS,1)';
DS=minObjV;
% 画出迭代图
GCLE=GCL(:,end);

function CH=Select(Chrom,ObjV,Gg)
% 选择操作,一般选择前20%的优秀个体进行保留
% 默认的优秀率
if nargin<=2.5 | isempty(Gg), Gg=0.2; end
% 种群数目
N=size(Chrom,1);
% 优秀个体数目
NS=max(floor(N*Gg+0.5),2);
% 优秀编号
[~,OI]=sortrows(ObjV,1);
po=OI(1:NS);
% 被选择的个体
CH=Chrom(po,:);
return

function CH=Reverse(Chrom,D)
% 逆转操作,所有个体均进行逆转操作,
% 种群数目N和基因长度L
[N,L]=size(Chrom);
% 计算路径长度
ObjV=PathLength(Chrom,D);
% 逆转操作
CH=Chrom;
for n=1:N
    % 需要逆转的一段字符的起止两个位置
    rn=randperm(L,2);
    r1=min(rn);
    r2=max(rn);
    CH(n,r1:r2)=CH(n,r2:-1:r1);
end
% 计算逆转之后的路径长度
OV=PathLength(CH,D); 
IN=(OV>=ObjV);
% 更换逆转之后的结果
CH(IN,:)=Chrom(IN,:);
return

function CH=Recombin(Chrom,D)
% 重组操作，所有个体均进行重组操作,
% 种群数目N和基因长度L
[N,L]=size(Chrom);
% 计算路径长度
ObjV=PathLength(Chrom,D);
% 重组操作
CH=Chrom;
for n=1:N
    % 需要重组的当前个体
    ch=Chrom(n,:);
    % 需要重组的一段字符的起止两个位置
    rn=randperm(L,2);
    r1=min(rn);
    r2=max(rn);
    % 重组操作
    chtemp=ch([r1:r2]);
    temp=chtemp(1);
    tp=temp;
    chre=chtemp([2:end]);
    chlen=length(chre);
    % 按照城市之间的距离远近进行操作
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
% 计算重组之后的路径长度
OV=PathLength(CH,D);
IN=(OV>=ObjV);
% 更换重组之后的结果
CH(IN,:)=Chrom(IN,:);
return

function Len=PathLength(Chrom,D)
% 计算环路长度
% 将城市的顺序补成环路
Chrom=[Chrom,Chrom(:,1)];
% 城市的数目
N=size(D,1);
% 种群的数目
M=size(Chrom,1);
% 为每个个体的长度开辟空间
Len=zeros(M,1);
% 计算长度
for m=1:M
    ch=Chrom(m,:);
    ch1=ch(1:end-1);
    ch2=ch(2:end);
    d=D(ch1,ch2);
    Len(m)=sum(diag(d));
end
return

function CH=Mutate(Chrom,Pm)
% 变异或换位操作，一般选择20%的个体进行操作
% 默认的变异概率
if nargin<=1.5 | isempty(Pm), Pm=0.2; end
% 种群数目N和基因长度L
[N,L]=size(Chrom);
% 为变异之后的新个体开辟存储空间
mn=round(Pm*N);         % 需要变异的个体数目
mC=randperm(N,mn);      % 变异的个体
CH=zeros(mn,L);         % 变异后的新个体
% 进行变异
for n=1:mn
    % 需要变异的当前个体
    mc=mC(n);
    ch=Chrom(mc,:);
    % 需要换位的两个字符的位置
    rn=randperm(L,2);
    r1=rn(1);
    r2=rn(2);
    % 具体换位的两个字符
    Sr1=ch(r1);
    Sr2=ch(r2);
    % 进行换位
    ch(r1)=Sr2;
    ch(r2)=Sr1;
    % 保存结果
    CH(n,:)=ch;
end
return

function CH=Gression(Chrom,D)
% 移位操作，所有个体均进行移位操作,
% 种群数目N和基因长度L
[N,L]=size(Chrom);
% 计算路径长度
ObjV=PathLength(Chrom,D);
% 进行移位
CH=Chrom;
for n=1:N
    % 需要移位的当前个体
    ch=Chrom(n,:);
    % 需要移位的一段字符的起止两个位置
    rn=randperm(L,2);
    r1=min(rn);
    r2=max(rn);
    % 字符串移位的操作
    chtemp=ch([(r1+1):r2]);
    chtemp=[chtemp,ch(r1)];
    CH(n,r1:r2)=chtemp;
end
% 计算重组之后的路径长度
OV=PathLength(CH,D);
IN=(OV>=ObjV);
% 更换重组之后的结果
CH(IN,:)=Chrom(IN,:);
return

function CH=GAcross(Chrom,Pc)
% 交叉操作的主程序，一半选择80%的个体组进行操作
% 默认的交叉概率
if nargin<=1.5 | isempty(Pc), Pc=0.8; end
% 种群数目
[N,L]=size(Chrom);
% 相互交叉的两个个体
[X,Y]=meshgrid([1:N]);
XY=[X(:),Y(:)];
% 删除两个个体相同的情况
po=find(abs(XY(:,1)-XY(:,2))<1e-8);
XY(po,:)=[];
% 可以交叉的总个体数目
XYN=size(XY,1);
% 为交叉之后的新个体开辟存储空间
mn=round(Pc*N);         % 需要交叉的个体数目
mC=randperm(XYN,mn);    % 交叉的个体组
CH=zeros(2*mn,L);       % 交叉之后的新个体
% 进行交叉
for n=1:mn
    % 需要交叉的当前个体组
    mc=mC(n);
    r1=XY(mc,1);
    r2=XY(mc,2);
    ch1=Chrom(r1,:);
    ch2=Chrom(r2,:);
    % 交叉之后的结果
    [A,B]=intercross(ch1,ch2);
    CH(2*n-1,:)=A;
    CH(2*n,:)=B;
end
return

function [an,bn]=intercross(a,b)
% 交叉操作的子程序
% 城市的数目
L=length(a);
% 交叉的两个位置
rn=randperm(L,2);
% 交叉的部分(从r1开始,到r2终止)
if abs(rn(1)-rn(2))>0.5
    r1=min(rn);
    r2=max(rn);
else
    r1=rn(1);
    r2=L;
end
r12=r2-r1+1;
% 交叉之前的保留
a0=a;
b0=b;
% 相互冲突的部分
ab0=intersect(a0(r1:r2),b0(r1:r2));
abua=setdiff(b0(r1:r2),ab0);
abub=setdiff(a0(r1:r2),ab0);
% 交叉之后的
an=a;
bn=b;
an(r1:r2)=zeros(r12,1);
bn(r1:r2)=zeros(r12,1);
% 寻找冲突的部分
for m=1:length(abua)
    po=find(abs(an-abua(m))<1e-8);
    an(po)=abub(m);
end
for m=1:length(abub)
    po=find(abs(bn-abub(m))<1e-8);
    bn(po)=abua(m);
end
% 生成交叉之后的
an(r1:r2)=b0(r1:r2);
bn(r1:r2)=a0(r1:r2);
return
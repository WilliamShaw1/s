function W=Adm2Wvm(A)
% 将邻接矩阵A转化为权值矩阵W(有向图或无向图均可)

% 节点数目
N=size(A,1);
W=A;
[m,n]=find(abs(W)<1e-8);
% 无边节点间距
for k=1:length(m)
    if abs(m(k)-n(k))>1e-8
        W(m(k),n(k))=Inf;
    end
end
% 同一节点间距
for k=1:N
    W(k,k)=0;
end
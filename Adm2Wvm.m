function W=Adm2Wvm(A)
% ���ڽӾ���Aת��ΪȨֵ����W(����ͼ������ͼ����)

% �ڵ���Ŀ
N=size(A,1);
W=A;
[m,n]=find(abs(W)<1e-8);
% �ޱ߽ڵ���
for k=1:length(m)
    if abs(m(k)-n(k))>1e-8
        W(m(k),n(k))=Inf;
    end
end
% ͬһ�ڵ���
for k=1:N
    W(k,k)=0;
end
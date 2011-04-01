function [edien]=bcedge(ien,bcn,elementype) % decide the local node index 
% for elements on solid boundary 

n=elementype;
edien=zeros(1,n);
for i=1:n
    node=ien(i);
    indx=find(node == bcn);
    m=length(indx);
    if (m ~= 0)
        edien(i)=1;
    end 
end

        




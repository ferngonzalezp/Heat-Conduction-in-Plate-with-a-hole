 tk = 79.5;
 l = 6;
 h = 6; 
 li = 2;
 hi = 2;
 t = 1;
 nel = 512;
 dofs = 1;
 
 t1 = thermal2d(tk,l,h,t,hi,li);
 
 Conect = t1.conect_sqr(nel);
 Nodes = Conect.Nodes;
 Elem = Conect.Elem;
 
 %condiciones de contorno
 r1 = find(Nodes(:,1) == li & Nodes(:,2) >= hi & Nodes(:,2) <= h-hi);
 r2 = find(Nodes(:,2) == hi & Nodes(:,1) > li & Nodes(:,1) <= l-li);
 r3 = find(Nodes(:,1) == l-li & Nodes(:,2) > hi & Nodes(:,2) <= h-hi);
 r4 = find(Nodes(:,2) == h-hi & Nodes(:,1) > li & Nodes(:,1) < l-li);
 r5 = find(Nodes(:,1) == 0);
 r6 = find(Nodes(:,2) == 0 & Nodes(:,1) > 0);
 r7 = find(Nodes(:,1) == l & Nodes(:,2)> 0 & Nodes(:,2) < h);
 r8 = find(Nodes(:,2) == h & Nodes(:,1)> 0 & Nodes(:,1) <= l);

 ri = [r1;r2;r3;r4];
 re = [r5;r6;r7;r8];
 desp = [ri zeros(size(ri)) 50.*ones(size(ri))];
 desp = [desp; re zeros(size(re)) 30.*ones(size(re))];
 
 K = t1.stiffness(nel,desp,dofs);
 spy(K)
 figure(2)
 plot(Nodes(:,1),Nodes(:,2),'o')

 T = sparse(size(Nodes,1),1);
 T(ri) = 50;
 T(re) = 30;
 r = 1:size(Nodes,1);
 r([ri' re']) = [];
 
 K_unknown = K;
 K_unknown(:,[ri;re]) = [];
 
T(r') = K_unknown\(-K(:,[ri;re])*T([ri;re]));
 
    figure(3)
for i = 1:size(Elem,1)
    patch(Nodes(Elem(i,:),1)',Nodes(Elem(i,:),2)',T(Elem(i,:))')
    hold on
end
    colorbar

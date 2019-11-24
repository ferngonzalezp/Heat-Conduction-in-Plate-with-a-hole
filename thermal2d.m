classdef thermal2d 
    properties
        tk
        l 
        h
        li
        hi
        t
    end
    methods
        function V = conect_sqr(obj,nel)
            nel1 = nel;
                while(1)
                    if mod(nel,2) == 0
                        nx = nel/2;
                    else
                        nx = nel;
                    end
                    ny = nel/2;
                    if 8*nx*ny == nel1
                        break;
                    else
                        nel = nx;
                    end
                end
                nel = 8*nx*ny;
                x = linspace(0, obj.l,3*nx+1);
                y = linspace(0,obj.h,3*ny+1);
                [X,Y] = meshgrid(x,y);
                Nodes = [X(:) Y(:)];
                Elem = zeros(nel,4);
                k = find(Nodes(:,1) > (obj.l-obj.li)/2 & Nodes(:,1) < (obj.l-obj.li)/2 + obj.li & Nodes(:,2) > (obj.h-obj.hi)/2 & Nodes(:,2) < (obj.h-obj.hi)/2 + obj.hi );
                Nodes(k,:) = [];
                node1 = 1;
            for i = 1:nel
                    Elem(i,1) = node1;
                    Elem(i,2) = find(Nodes(:,1) == Nodes(Elem(i,1),1) & Nodes(:,2) > Nodes(Elem(i,1),2),1);
                    Elem(i,4) = find(Nodes(:,2) == Nodes(Elem(i,1),2) & Nodes(:,1) > Nodes(Elem(i,1),1),1);
                    Elem(i,3) = find(Nodes(:,1) == Nodes(Elem(i,4),1) & Nodes(:,2) > Nodes(Elem(i,4),2),1);
                    if Nodes(Elem(i,3),2) == obj.l
                     node1 = find(Nodes(:,1) > Nodes(Elem(i,1),1),1);
                    else
                     if Nodes(Elem(i,2),1) >= obj.li && Nodes(Elem(i,2),2) >= obj.hi && Nodes(Elem(i,2),1) < (obj.l-obj.li) && Nodes(Elem(i,2),2) < (obj.h-obj.hi)
                        node1 = find(Nodes(:,1) == Nodes(Elem(i,1),1) & Nodes(:,2) == (obj.h-obj.hi)/2+obj.hi,1);
                     else
                         node1 = Elem(i,2);
                     end
                    end
            end
            V.Nodes = Nodes;
            V.Elem = Elem;
        end
        function V = dof_list(obj,dofs,desp,nel)
            Nodes = obj.conect_sqr(nel).Nodes;
            dof_fixed = size(Nodes,1)*dofs*10;
            dof_free = 1;
            V.total_dof = size(Nodes,1)*dofs*10;
            V.dof_list = zeros(size(Nodes,1),dofs+1);
            V.dof_list(:,1)=(1:1:size(Nodes,1));
            for i=1:size(desp,1)
            if desp(i,2) == 1
                a = dof_fixed;
                dof_fixed = dof_fixed + 1;
            else
                a = 0;
            end
            V.dof_list(desp(i,1),:) = [desp(i,1), a];
            end
        for i = 1:size(V.dof_list,1)

                if V.dof_list(i,2) == 0
                    a = dof_free;
                    dof_free = dof_free+1;
                    V.dof_list(i,2) = a;
                end
        end
        V.dof_free = dof_free;
        end
        function v = material(obj)
            D = [obj.tk 0; 0 obj.tk];
            v = D;
        end
        function v = shapefun1(obj,chi,ep)
            v.f = [1/4*(1-chi)*(1-ep) 1/4*(1-chi)*(1+ep) 1/4*(1+chi)*(1+ep) 1/4*(1+chi)*(1-ep)];
            v.df = [  ep/4 - 1/4,    -1/4 - ep/4,  ep/4 + 1/4, -ep/4 + 1/4; chi/4 - 1/4, - chi/4 + 1/4, chi/4 + 1/4,  -1/4 - chi/4];
        end
        function v = integrate2(obj,f)
            v = f(-sqrt(1/3),-sqrt(1/3))+f(-sqrt(1/3),sqrt(1/3))+f(sqrt(1/3),-sqrt(1/3))+f(sqrt(1/3),sqrt(1/3));
        end
        function v = strain(obj,chi,ep,i,nel)
            Nodes = obj.conect_sqr(nel).Nodes;
            Elem = obj.conect_sqr(nel).Elem;
            N = obj.shapefun1(chi,ep);
            Je = [0 0; 0 0];
            for j = 1:length(N.f)
                Je = Je + [N.df(1,j)*Nodes(Elem(i,j),1) N.df(1,j)*Nodes(Elem(i,j),2); N.df(2,j)*Nodes(Elem(i,j),1) N.df(2,j)*Nodes(Elem(i,j),2)];
            end
            Be = [];
            for j = 1:length(N.f)
                Be = [Be  Je\N.df(:,j)];
            end
            v.Be = Be;
            v.Je = Je;
        end
        function v = stiffness(obj,nel,desp,dofs)
            dof = obj.dof_list(dofs,desp,nel);
            dof_free = dof.dof_free;
            D = obj.material();
            K = sparse(dof_free-1, dof_free-1);
            Elem = obj.conect_sqr(nel).Elem;
            for i = 1:size(Elem,1)
                conect = index(Elem(i,:),dofs,dof.dof_list);
                f = @(chi,ep) obj.strain(chi,ep,i,nel).Be'*D*obj.strain(chi,ep,i,nel).Be*obj.t*det(obj.strain(chi,ep,i,nel).Je);
                    Ke = obj.integrate2(f);
                for k = 1: size(Ke,1)
                                if conect(k) < dof.total_dof
                                    for j = 1: size(Ke,2)
                                        if conect(j) < dof.total_dof
                                        K(conect(k),conect(j)) = Ke(k,j)+ K(conect(k),conect(j));
                                        end
                                    end
                                end
                end
            end
            v = K;
        end
        function obj = thermal2d(tk,l,h,t,hi,li)
            obj.tk = tk;
            obj.l = l;
            obj.h = h;
            obj.t = t;
            obj.hi = hi;
            obj.li = li;
        end
    end
end
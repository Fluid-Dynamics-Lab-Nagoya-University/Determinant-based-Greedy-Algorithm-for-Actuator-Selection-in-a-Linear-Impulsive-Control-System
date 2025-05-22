%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%F_act_dg_under1.m ver.240621%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [obj_select,index_select]=F_act_dg_under1(W,p)
    r=size(W,2);
    index_select=zeros(p,1);
    obj_select=zeros(p,1);
    Ws=zeros(p,r);
    WWAI=zeros(p,p);
    WWAI_new=zeros(p,p);
    for k=1:p
        if k==1
            obj=sum(abs(W).^2,2);
        else
            F=eye(r)-Ws(1:k-1,:)'*WWAI(1:k-1,1:k-1)*Ws(1:k-1,:);
            obj=sum((W*F).*conj(W),2);
            obj=abs(obj);
            obj(index_select(1:k-1,1),1)=-Inf;
        end
        [obj_max,index_obj_max]=max(obj);
        index_select(k,1)=index_obj_max;
        obj_select(k,1)=obj_max;
        Ws(k,:)=W(index_obj_max,:);
        if k==1
            WWAI(1,1)=1/(Ws(1,:)*Ws(1,:)');
        else
            denominator=Ws(k,:)*(eye(r)-Ws(1:k-1,:)'*WWAI(1:k-1,1:k-1)*Ws(1:k-1,:))*Ws(k,:)';
            WWAI_new(1:k-1,1:k-1)=WWAI(1:k-1,1:k-1)*(eye(k-1)+Ws(1:k-1,:)*Ws(k,:)'*Ws(k,:)*Ws(1:k-1,:)'*WWAI(1:k-1,1:k-1)/denominator);
            WWAI_new(k,1:k-1)=-Ws(k,:)*Ws(1:k-1,:)'*WWAI(1:k-1,1:k-1)/denominator;
            WWAI_new(1:k-1,k)=WWAI_new(k,1:k-1)';
            WWAI_new(k,k)=1/denominator;
            WWAI=WWAI_new;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
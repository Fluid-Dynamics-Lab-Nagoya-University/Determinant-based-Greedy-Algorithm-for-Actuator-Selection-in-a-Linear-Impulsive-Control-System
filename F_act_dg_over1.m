%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%F_act_dg_over1.m ver.240621%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [obj_select,index_select]=F_act_dg_over1(W,p,obj_select_under,index_select_under)
    r=size(W,2);
    index_select=zeros(p,1);
    index_select(1:r,1)=index_select_under;
    obj_select=zeros(p,1);
    obj_select(1:r,1)=obj_select_under;
    Ws=zeros(p,r);
    Ws(1:r,:)=W(index_select_under,:);
    WAWI=inv(Ws(1:r,:)'*Ws(1:r,:));
    for k=r+1:p
        obj=sum((W*WAWI).*conj(W),2);
        obj=abs(obj);
        obj(index_select(1:k-1,1),1)=-Inf;
        [obj_max,index_obj_max]=max(obj);
        index_select(k,1)=index_obj_max;
        obj_select(k,1)=obj_max;
        Ws(k,:)=W(index_obj_max,:);
        denominator=1+Ws(k,:)*WAWI*Ws(k,:)';
        WAWI=WAWI*(eye(r)-Ws(k,:)'*Ws(k,:)*WAWI/denominator);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfrep=I_max_new;
%dis=1;
im_bin1 = zeros(size(tfrep));
a=size(im_bin1);
T=0.02*max(max(tfrep));
for kk = 1:a(2)
    q1 = tfrep(:,kk);
    u = diff(q1);
    uu = q1;
    zc_max =[];
    count_max = 1;
    for zz = 1:(length(u)-2)
        if u(zz)>=0
            if u(zz+1)<0
                if q1(zz)>T
                    zc_max(count_max) = zz;
                    count_max = count_max+1;
                end
            end
        end
    end
    q2 = zeros(1,a(1));
    q2(zc_max+1) = 1;
    im_bin1(kk,:) = q2;
    clear q* z* c* ref
end
M=10;
ref = [0:M+1 -1:-1:-(M+1)];
[val, idx] = sort(abs(ref));
nref = ref(idx);
sch = [ones(1,length(nref)); nref];
[el,ei] = edge_link(im_bin1, length(x)*2/3, sch);


% if length(el)==1
%     %mse(2)=sum(IF_O(:,1).^2);
%     mse(2)=1;
%     mse(3)=1;
% elseif  isempty((el))
%     %mse(1)=sum(IF_O(:,2).^2);
%     %mse(2)=sum(IF_O(:,1).^2);
%     mse(1)=1;
%     mse(2)=1;
%     mse(3)=1;
% elseif length(el)==2
%     mse(3)=1;
% 
% end

% if length(el)>N_C
%     mse=[];
%     
%     for i=1:N_C
%         [ value index]=min(msee);
%         msee(index(1))=1000;
%         mse(i)=value;
%     end
% elseif length(el)==N_C
%     mse=msee;
% end


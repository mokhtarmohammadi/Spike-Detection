tfrep=I9;
%dis=1;
im_bin1 = zeros(size(tfrep));
aa=size(im_bin1);
T=0.1*max(max(tfrep));

for kk = 1:aa(2)
    q1 = tfrep(:,kk);
    u = diff(q1);
    uu = q1;
    zc_max =[];
    count_max = 1;
    for zz = 1:(length(u)-2)
        if u(zz)>=0
            if u(zz+1)<0
                if q1(zz)>T*max(q1)
                    zc_max(count_max) = zz;
                    count_max = count_max+1;
                end
            end
        end
    end
    q2 = zeros(1,aa(1));
    q2(zc_max+1) = 1;
    im_bin1(kk,:) = q2;
    clear q* z* c* ref
end
MM=10;
ref = [0:MM+1 -1:-1:-(MM+1)];
[val, idx] = sort(abs(ref));
nref = ref(idx);
sch = [ones(1,length(nref)); nref];
[el,ei] = edge_link(im_bin1,2*length(s)/3, sch);
dd=ei';


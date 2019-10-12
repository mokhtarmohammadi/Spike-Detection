dis=0;


LL=22;
IF1=zeros([N_S length(s) N_C]);

%  for k1=1:N_S
slope=0;
%I_max=I_max_new(end:-1:1,:);
I_max=I;
I_max2=I_max;
I_max1=I_max;
IF=zeros(length(s),M);
%LL=9;
for k=1:N_C
    I_max=I_max1;
    [xx yy]=find(I_max1(:,30:100)==max(max(I_max1(:,30:100))));
    if and((xx(1)-LL)>0,xx(1)+LL<length(s))
        I_max1(xx(1)-LL:xx(1)+LL,yy)=0;
    else
        I_max1(xx(1),yy)=0;
    end
    
    IF(yy,k)=xx(1);
    x=xx;
    x_prev=x;
    for i=yy+1:length(s)
        %l=zeros(1,length(xn));
        slope=0.1*(x-x_prev)+0.9*slope;
        slope=0;
        x_prev=x;
        x=x+slope;
        x=round(x);
        
        if and((x(1)-LL)>0,x(1)+LL<length(s))
            I_max(x(1)-LL:x(1)+LL,i)=100*I_max(x(1)-LL:x(1)+LL,i);
        elseif and((x(1))>0,x(1)<length(s))
            
            I_max(x(1),i)=100*I_max(x(1),i);
        else
            x=x-slope;
            x=round(x);
            I_max(x(1),i)=100*I_max(x(1),i);
        end
        
        
        
        %x_prev=x;
        [x y]=find(I_max(:,i)==max(I_max(:,i)));
        x=mean(x);
        if and((x(1)-LL)>0,x(1)+LL<length(s))
            I_max1(x(1)-LL:x(1)+LL,i)=0;
        else
            I_max1(x(1),i)=0;
        end
        
        
        
        
        
        IF(i,k)=x(1);
    end
    
    
    x=xx;
    x_prev=x;
    for i=yy-1:-1:1
        %l=zeros(1,length(s));
        slope=(x-x_prev)*0.1+0.9*slope;
        slope=0;
        x_prev=x;
        x=x+slope;
        x=round(x);
        if and((x(1)-LL)>0,x(1)+LL<length(s))
            I_max(x(1)-LL:x(1)+LL,i)=100*I_max(x(1)-LL:x(1)+LL,i);
        else
            I_max(x(1),i)=100*I_max(x(1),i);
        end
        
        
        
        %x_prev=x;
        
        [x y]=find(I_max(:,i)==max(I_max(:,i)));
        x=mean(x);
        if and((x(1)-LL)>0,x(1)+LL<length(s))
            I_max1(x(1)-LL:x(1)+LL,i)=0;
        elseif and((x(1))>0,x(1)<length(s))
            
            I_max(x(1),i)=100*I_max(x(1),i);
        else
            x=x-slope;
            I_max(x(1),i)=100*I_max(x(1),i);
        end
        
        
        
        
        
        IF(i,k)=x(1);
    end
    
    if dis==1
      figure;
                 imshow(1-I_max1,[])
    end
end
IFF=IF;


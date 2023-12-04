function[h]=histogram(datos)
datos=datos(:);
ind=find(isnan(datos)==1);
datos(ind)=0;
ind=find(isinf(datos)==1);
datos(ind)=0;
tam=length(datos);
m=ceil(max(datos))+1;
h=zeros(1,m);
for i=1:tam,
    f=floor(datos(i));    
    if(f>0 & f<(m-1))        
        a2=datos(i)-f;
        a1=1-a2;
        h(f)  =h(f)  + a1;      
        h(f+1)=h(f+1)+ a2;                          
    end;
end;
h=conv(h,[1,2,3,2,1]);
h=h(3:(length(h)-2));
h=h/sum(h);

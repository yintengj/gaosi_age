%% ��˹���ģ�����㷨
%% ���룺�������ݼ�x��ͳ�Ʒֲ�Ƶ��h�����ȡֵm����Сȡֵmi����˹�ֲ�����k
%% �������˹�ֲ���������׼���ֵ�����ʡ�����
function [mu,v,p,prb,k]=GMM(x,h,m,mi,k)

 if (k == 0)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �û���ָ������������������ֵ����Ϊ��ʼ��˹�ֲ��������˲��ֳ���ִ���ɣ������и��õķ���
    for i = 2:(size(x,1)-1)
        if (h(i)>h(i-1))&&(h(i)>h(i+1))
            k = k+1;
            mu(k) = x(i);
        end
    end
%% ��ֵ����̫�࣬һ�ִ�����ȥ����������
    l = 1
    n = k
    while(l<n)
        if((mu(l+1)-mu(l))<(m-mi)/(k+2))
            mu(l+1)=[]
            n = n-1
        else
            l = l+1;
        end
    end
    k = l;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mu = (1:k)*m/(k+1);         % ����ֲ���С����ֵ���������
    v=ones(1,k)*m;
 else 
    mu = (1:k)*(m-mi)/(k+1)+mi; % ����ֲ�û�й��ȼ��и��ݳ�ʶ�ǽ��������ݽ���k�ȷ�
    v=ones(1,k)*(m-mi);
 end
 p=ones(1,k)*1/k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������������Ͽ�Դ��
sml = mean(diff(x))/1000;
nIter = 0;
while(1)
        % Expectation
        nIter = nIter+1;
        prb = distribution(mu,v,p,x);
        scal = sum(prb,2)+eps;
        loglik=sum(h.*log(scal));
        
        %Maximizarion
        for j=1:k
                pp=h.*prb(:,j)./scal;
                p(j) = sum(pp);
                mu(j) = sum(x.*pp)/p(j);
                vr = (x-mu(j));
                v(j)=sum(vr.*vr.*pp)/p(j)+sml;
        end
        p = p + 1e-3;
        p = p/sum(p);
        % Exit condition
        prb = distribution(mu,v,p,x);
        scal = sum(prb,2)+eps;
        nloglik=sum(h.*log(scal)); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp = 0.001;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if((nloglik-loglik)<temp) 
            break;
        end;        

        figure(1)
        clf
        plot(x,h);
        hold on   
        plot(x,prb,'g--')
        temp = sum(prb,2);
        plot(x,temp,'r')
%       plot(x,abs(h-temp),'k')
        XLabel('');
        YLabel('');
        drawnow
end
v = sqrt(v);








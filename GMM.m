%% 高斯混合模型主算法
%% 输入：样本数据集x、统计分布频度h、最大取值m、最小取值mi、高斯分布数量k
%% 输出：高斯分布期望、标准差、幅值、概率、数量
function [mu,v,p,prb,k]=GMM(x,h,m,mi,k)

 if (k == 0)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 用户不指定概念数量，则计算峰值点作为初始高斯分布数量，此部分程序仓促完成，可能有更好的方法
    for i = 2:(size(x,1)-1)
        if (h(i)>h(i-1))&&(h(i)>h(i+1))
            k = k+1;
            mu(k) = x(i);
        end
    end
%% 峰值数量太多，一种粗蛮的去掉抖动方法
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
    mu = (1:k)*m/(k+1);         % 如果分布在小的数值区域更集中
    v=ones(1,k)*m;
 else 
    mu = (1:k)*(m-mi)/(k+1)+mi; % 如果分布没有过度集中根据常识是将论域数据进行k等分
    v=ones(1,k)*(m-mi);
 end
 p=ones(1,k)*1/k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 下面这段是网上开源的
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








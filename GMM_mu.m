%% 指定高斯分布数量和期望的高斯混合模型主算法
%% 输入：样本数据集x、统计分布频度h、最大取值m、最小取值mi、高斯分布数量k、高斯分布期望
%% 输出：高斯分布期望、标准差、幅值、概率

function [mu,v,p,prb]=GMM_mu(x,h,m,mi,k,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 开源代码
 v=ones(1,k)*m;
 p=ones(1,k)*1/k;
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
        temp = sum(prb,2)
        plot(x,temp,'r')
        plot(x,abs(h-temp),'k')
        XLabel('');
        YLabel('');
        drawnow
end
 v = sqrt(v);








%% 从高斯分布的交叠程度计算高斯云的熵、超熵和概念含混度
%% 输入：高斯分布的期望数组、标准差数组、幅值数组
%% 输出：高斯云的熵数组、超熵数组、概念含混度数组
function [En,He,belta] = comp_En_He(mu,v,p)

mu=mu(:);
v=v(:);
p=p(:);
k = size(mu,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 两个以上高斯分布才可能存在交叠现象
if k>1     
    for i = 1:k-1
       if ((mu(i)+3*v(i))<=(mu(i+1)-3*v(i+1)))% 如果弱外围元素不相交
           temp_alpha(i) = 1;
           En(i) = v(i);
           He(i) = 0;    % 高斯变换生成的高斯分布的外围节点与其他高斯分布不相交可以独立表示一个概念，可以忽略粒度的不确定性。
       else
           temp_alpha(i) = [mu(i+1)-mu(i)]/[3*v(i+1)+3*v(i)];% 保证弱外围元素不相交的边界条件
           He(i) = (1-temp_alpha(i))*v(i)/6 ;% 方差变化范围[temp*v,v],从弱外围元素不相交到骨干元素不相交 
           En(i) = v(i)-3*He(i); 
       end
       if i>1    
            alpha(i) = min(temp_alpha(i),temp_alpha(i-1));
            He(i) = (1-alpha(i))*v(i)/6 ;% 方差变化范围[temp*v,v],从弱外围元素不相交到骨干元素不相交 
            En(i) = v(i)-3*He(i);      
       end
    end     
       if (He(k-1) ~= 0) % 判断第如果弱外围元素相交
            He(k) = (1-temp_alpha(k-1))*v(k)/6 ;
            En(k) = v(k)-3*He(k);
       else
           En(k) = v(k);
           He(k) = 0;
       end
elseif k ==1
     En(k) = v(k);
     He(k) = v(k)/3;
end
for l=1:k
   belta(l) = 3*He(l)/En(l);
end
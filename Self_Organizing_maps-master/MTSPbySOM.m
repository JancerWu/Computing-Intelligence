clc
clear all
close all
%% 输入文件
%tspdata=importdata('att48.tsp.txt')
%tspdata=importdata('bayg29.tsp.txt')
%tspdata=importdata('eil76.tsp.txt')
tspdata=importdata('eil101.tsp.txt');
% tspdata=importdata('pa561.tsp.txt');
%tspdata=importdata('tsp225.tsp.txt');
%% 获取图像的中心位置，将图像从中心位置建立XOY坐标系
datanew = tspdata(:, 2:3); % 获取训练集中的横纵坐标，根据横纵坐标获取中心点
maxv = max(datanew); % 获取横坐标的最大值，和纵坐标的最大值
minv = min(datanew);
maxvalue = maxv(1)*maxv(1)+maxv(2)*maxv(2);
wcenter = (maxv - minv)/2 + minv; % 初始权值的中心,中心有横纵坐标
%% 参数设定
alpha = 0.05; % 该学习率对结果的影响，越大训练速度越快，但准确率越低 [0.001，0.1]
beta = 0.1; % 学习率 [0.001，0.5]
gain = 10; % 可以是任意的，10、15、20
percent = 0.5; % 邻居的占比
saleman = 4; %旅行商数量
nsize = size(tspdata); % size返回行数和列数
ncity = nsize(1); % 城市个数，即nsize的第1列
m = ncity;
k_bias = 0.75;
winit = [rands(m, 2)+wcenter+wcenter.*[0,k_bias],rands(m, 2)+wcenter+wcenter.*[0,-k_bias]...
    ,rands(m, 2)+wcenter+wcenter.*[k_bias,0],rands(m, 2)+wcenter+wcenter.*[-k_bias,0]];
w = winit;
wold = w;  %用来判断权值是否和上一个相等以便退出循环
times = 0; %迭代次数
inhibit = zeros(m,1); % 所有结点的抵达状态, m代表城市个数
min_x = 0;
min_y = 0;
%% mindst寻找图像中最近两点的距离，为了使函数在各个尺度的数据中都用较高的效率
[mindst] = findmin_distance(datanew);

while times<150
%% 移动点      
        
    plot(datanew(:,1), datanew(:,2), 'ko', 'MarkerFaceColor', 'r');
    title('solution based on som algorithm');
    hold on;
    plot(w(:,1), w(:,2), '.-');
    plot(w(:,3), w(:,4), '.-');
    plot(w(:,5), w(:,6), '.-');
    plot(w(:,7), w(:,8), '.-');
%     for i=1:saleman
%         plot(w(:,1+(2*(salemannum-1))), w(:,2+(2*(salemannum-1))), '.-');
%     end
    pause(0.0001); % 暂停0.0001秒后执行下一条指令
    hold off;
    inhibit = zeros(m,1); % 重置所有结点的抵制状态, m代表城市个数
    %随机输入城市的编号
    yrand = randperm(ncity);
    for choice_num = 1:ncity
        newidx = yrand(choice_num);
        a = datanew(newidx, :); % 获取城市的坐标
        % 计算结点(1..m) 和 各个网络间的距离
        for salemannum=1:saleman
            for j=1:m
                if inhibit(j) == 1
                    T(salemannum,j) = maxvalue;%如果这个点已经能找到最近点，就不再计算它
                else
                    T(salemannum,j) = (a(1)-w(j,1+(2*(salemannum-1))))^2 + (a(2)-w(j,2+(2*(salemannum-1))))^2;
                end
            end
        end
        % [Tmin保存这个获胜节点到当前城市的短距离， min_x表示是离那个网络最近，min_y表示那个城市的下标]
        Tmin = min(min(T));
        [min_x,min_y]=find(Tmin == T);

        inhibit(min_y) = 1;
        % 将每次获胜的节点保存下来
        y(min_y) = newidx;
        f = zeros(1,m);
        %只更新离当前网络最近的点所在的网络
         for j=1:m
            % 表示在环上节点 j 和 J 的距离
            d=min(abs(j-min_y), m-abs(j-min_y));
            [size_dx, size_dy]= size(d);
            if size_dx>=2
                d = d(1,1);
            end
            if d<percent*m
                f(j) = exp(-d*d/(gain*gain));
                w(j,1+(2*(min_x-1)))=w(j,1+(2*(min_x-1)))+beta*f(j)*(a(1)-w(j,1+(2*(min_x-1))));
                w(j,2+(2*(min_x-1)))=w(j,2+(2*(min_x-1)))+beta*f(j)*(a(2)-w(j,2+(2*(min_x-1))));
            end
        end
        
        % 保存竞争胜利的结点如果距离小于0.01*mindst
        distJ = sqrt((a(1)-w(min_y,1+(2*(salemannum-1)))).^2 + (a(2)-w(min_y,2+(2*(salemannum-1)))).^2);
        if distJ < 0.01*mindst
            w(min_y,1+(2*(salemannum-1)))=a(1);
            w(min_y,2+(2*(salemannum-1)))=a(2);
        end
    end
    % 学习速率改变
    alpha = alpha*0.998;
    gain = (1-alpha)*gain;
    
    % 程序退出条件，当权值不再改变时
     if w==wold
        break;
     end
    wold = w;
    times = times + 1;
    
end

%delete
k = 1;
z = zeros(1,ncity);
l =zeros(1,saleman);
for salemannum=1:saleman
    for i=1:m
        for j = 1:ncity
            if(abs(norm(w(i,1+2*(salemannum-1):2*(salemannum))-datanew(j,:)))<2 && z(1,j)~=1)
                new_w(k,1+2*(salemannum-1):2*(salemannum)) =w(i,1+2*(salemannum-1):2*(salemannum));
                % 将每次获胜的节点保存下来
                y(k) = j;
                k = k +1;
                l(1,salemannum) = l(1,salemannum)+1;
                z(j)=1;
                break;
            end
        end
    end
end
l(2) = l(1)+l(2);
l(3) = l(3)+l(2);
l(4) = l(4)+l(3);
%% 画出路径 可以简化
plot(new_w(1:l(1,1),1), new_w(1:l(1,1),2), 'ko', 'MarkerFaceColor', 'b');
hold on;
plot(new_w(l(1,1)+1:l(2),3), new_w(l(1,1)+1:l(2),4), 'ko', 'MarkerFaceColor', 'r');
hold on;
plot(new_w(l(1,2)+1:l(3),5), new_w(l(1,2)+1:l(3),6), 'ko', 'MarkerFaceColor', 'g');
hold on;
plot(new_w(l(1,3)+1:l(4),7), new_w(l(1,3)+1:l(4),8), 'ko', 'MarkerFaceColor', 'y');
hold on;
title('TSP route by som ');
hold on;
plot(new_w(1:l(1,1),1), new_w(1:l(1,1),2), 'b.-');
plot([new_w(1,1) new_w(l(1,1),1) ],[new_w(1,2) new_w(l(1,1),2)], 'b.-');
hold on;
plot(new_w(l(1,1)+1:l(2),3), new_w(l(1,1)+1:l(2),4), 'r.-');
plot([new_w(l(1,1)+1,3) new_w(l(2),3)],[new_w(l(1,1)+1,4) new_w(l(2),4)], 'r.-');
hold on;
plot(new_w(l(1,2)+1:l(3),5), new_w(l(1,2)+1:l(3),6), 'g.-');
plot([new_w(l(1,2)+1,5) new_w(l(3),5)],[new_w(l(1,2)+1,6) new_w(l(3),6)], 'g.-');
hold on;
plot(new_w(l(1,3)+1:l(4),7), new_w(l(1,3)+1:l(4),8), 'y.-');
plot([new_w(l(1,3)+1,7) new_w(l(4),7)],[new_w(l(1,3)+1,8) new_w(l(4),8)], 'y.-');
hold on;
pause(0.0001); % 暂停0.0001秒后执行下一条指令
hold off;
% plot(new_w(:,3), new_w(:,4), 'b.-');
% hold on;
% plot(new_w(:,5), new_w(:,6), 'b.-');
% plot(new_w(:,7), new_w(:,8), 'b.-');
pause(0.0001); % 暂停0.0001秒后执行下一条指令
hold off;
% 输出迭代次数
times
% 保存路径
solution = y;

clc
clear all
close all
tic
%% 输入文件
% tspdata=importdata('att48.tsp.txt');
% tspdata=importdata('eil51.tsp.txt');
%tspdata=importdata('bayg29.tsp.txt');
% tspdata=importdata('eil76.tsp.txt');
tspdata=importdata('eil101.tsp.txt');
% tspdata=importdata('pa561.tsp.txt');
% tspdata=importdata('tsp225.tsp.txt');
%datanew = tspdata(:, 2:3); % 获取训练集中的横纵坐标和纵坐标
%% 获取图像的中心位置，将图像从中心位置建立XOY坐标系
datanew = tspdata(:, 2:3); % 获取训练集中的横纵坐标，根据横纵坐标获取中心点
maxv = max(datanew); % 获取横坐标的最大值，和纵坐标的最大值
minv = min(datanew);
maxvalue = maxv(1)*maxv(1)+maxv(2)*maxv(2);
wcenter = (maxv - minv)/2 + minv; % 初始权值的中心,中心有横纵坐标
%% 参数设定
alpha = 0.05; % 该学习率对结果的影响，越大训练速度越快，但准确率越低 [0.001，0.1]
beta = 0.2; % 学习率 [0.001，0.5]
gain = 10; % 可以是任意的，10、15、20
percent = 0.8; % 邻居的占比
saleman = 4; %旅行商数量
nsize = size(tspdata); % size返回行数和列数
ncity = nsize(1); % 城市个数，即nsize的第1列
m = ncity*2;
% m = ncity;
winit = [rands(m, 2)+wcenter];
w = winit;
wold = w;  %用来判断权值是否和上一个相等以便退出循环
times = 0; %迭代次数
new_w = ones(ncity, 2);%存放删除节点后的值
%% mindst寻找图像中最近两点的距离，为了使函数在各个尺度的数据中都用较高的效率
[mindst] = findmin_distance(datanew);

%%
while 1>0
    plot(datanew(:,1), datanew(:,2), 'ko', 'MarkerFaceColor', 'r');
    title('TSP route by som');
    hold on;
    plot(w(:,1), w(:,2), 'b.-');
    plot([w(1,1) w(m,1)], [w(1,2) w(m,2)], 'b.-');
    pause(0.0001); % 暂停0.0001秒后执行下一条指令
    hold off;
    inhibit = zeros(m,1); % 重置所有结点的抵制状态, m代表城市个数
    %随机输入城市的编号
    
    yrand = randperm(ncity);
    for pattern = 1:ncity
        newidx = yrand(pattern);
        a = datanew(newidx, :); % 获取城市的坐标
        % 计算结点(1..m) 和 输入结点a 的距离
        for j=1:m
            if inhibit(j) == 1
                T(j) = maxvalue;%如果这个点已经能找到最近点，就不再计算它
            else
                T(j) = (a(1)-w(j,1))^2 + (a(2)-w(j,2))^2;
            end
        end
        % [Tmin保存这个获胜节点到当前城市的短距离， Jmin保存这个获胜节点的下标]
        [Tmin, Jmin] = min(T);
        inhibit(Jmin) = 1;
        
        f = zeros(1,m);
        for j=1:m
            % 表示在环上节点 j 和 J 的距离
            d=min(abs(j-Jmin), m-abs(j-Jmin));
            % 对领域内的所有城市进行权值修改
            if d<percent*m
                f(j) = exp(-d*d/(gain*gain));
                w(j,1)=w(j,1)+beta*f(j)*(a(1)-w(j,1));
                w(j,2)=w(j,2)+beta*f(j)*(a(2)-w(j,2));
            end
        end
        
        % 保存竞争胜利的结点如果距离小于0.01*mindst
        distJ = sqrt((a(1)-w(Jmin,1))^2 + (a(2)-w(Jmin,2))^2);
        if distJ < 0.01*mindst
            w(Jmin,1)=a(1);
            w(Jmin,2)=a(2);
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
z = zeros(ncity);
for i=1:m
    for j = 1:ncity
        if(abs(norm(w(i,:)-datanew(j,:)))<0.1 && z(j)~=1)
            new_w(k,:) = w(i,:);
            % 将每次获胜的节点保存下来
            y(k) = j;
            k = k +1;
            z(j)=1;
            break;
        end
    end
end
plot(new_w(:,1), new_w(:,2), 'ko', 'MarkerFaceColor', 'r');
title('TSP route by som ');
hold on;
plot(new_w(:,1), new_w(:,2), 'b.-');
plot([new_w(1,1) new_w(ncity,1)], [new_w(1,2) new_w(ncity,2)], 'b.-');
pause(0.0001); % 暂停0.0001秒后执行下一条指令
hold off;

% 输出迭代次数
times
% 保存路径
solution = y;
solution
%计算旅行距离
tourdistance = 0;
for i=1:ncity-1
    tourdistance = tourdistance + norm(new_w(i+1)-new_w(i));
end
tourdistance = tourdistance + norm(new_w(ncity)-new_w(1));
tourdistance
toc
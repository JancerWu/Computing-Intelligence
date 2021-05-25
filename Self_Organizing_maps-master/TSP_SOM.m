clc
clear all
close all

% 导入城市数据,任选一
% cityData = load('eil101.tsp.txt');
% cityData = load('att48.tsp.txt');
% cityData = load('bayg29.tsp.txt');
% cityData = load('eil51.tsp.txt');
% cityData = load('eil76.tsp.txt');
cityData = load('pa561.tsp.txt');
% cityData = load('tsp225.tsp.txt');

% 获取城市坐标
city = cityData(: , 2 : 3);
% 计算城市样本中心

% 扩展为MTSP问题，根据环路数量计算多个聚类中心点
tourNum = 30;               % 环路数量
[~ , cityCenter] = kmeans(city, tourNum);
% 初始化参数
beta = 0.15;                % 更新权值的学习率1
neigborPercent = 0.6;      % 获胜节点的邻居占比
alpha = 0.04;              % 更新邻域的学习率2
g = 10;                    % 更新邻域的参数
cityNum = size(city, 1);
outNum = cityNum;
outW = zeros(tourNum, outNum, 2);  %最高维表示第几个环路
%初始化不同环路作图颜色分配方案，三原色随机取值
colorPack = rand(tourNum, 3);
for i = 1 : tourNum
    outW(i,:,:) = cityCenter(i) + rands(outNum, 2);
end
% outW = cityCenter + rands(cityNum, 2);    % 初始化权值
lastW = outW;   % 记录上一次的权值用于比较是否改变

% 以城市样本中两城市的最近距离作为基准，
% 后面训练中若权值距离小于该值的0.01，
% 则权值调整为该城市位置
mindst = inf;
for i = 1 : cityNum
    for j = i + 1 : cityNum
        dist = sqrt((city(i, 1) - city(j, 1)).^2 + (city(i, 2) - city(j, 2)).^2 );
        if (dist < mindst)
            mindst = dist;
        end
    end
end

% 开始训练
MaxG = 1000;   %最大迭代次数

% 计时开始
tic
stopCity = 0;
% cityOver = zeros(cityNum,1);
countTourCityNum = ones(tourNum,1);   %记录每个环路获胜节点个数
preserveNode = zeros(tourNum,cityNum,2);  %保留环路中的获胜节点
for time = 1 : MaxG
    
%   打印城市位置以及输出节点位置
    plot(city(:,1), city(:,2), 'or');
    hold on;
    for i = 1 : tourNum
        plot(outW(i,:,1), outW(i,:,2), '.-', 'Color', colorPack(i,:));
        plot([outW(i,1,1),outW(i,outNum,1)], [outW(i,1,2),outW(i,outNum,2)], '.-', 'Color', colorPack(i,:));
    end
    pause(0.001);
    hold off;
   
    done = zeros(cityNum, 1);  %更新所有输出节点的状态为初始状态
    numRandom = randperm(cityNum);   %随机给定一个城市顺序
    
%     recNum = 1;  %权值接近城市的节点数量
%     resTour = zeros(cityNum, 1); %近似解路径
%    计算给定城市位置与各权值的距离
    for oneByone = 1 : cityNum
        cityPos = city( numRandom(oneByone), :);
        dis = zeros(tourNum, outNum);
        for j = 1 : tourNum
            for i = 1 : outNum
                if done(i) == 1
                    dis(j,i) = inf;   
                else
                    dis(j,i) = (cityPos(1) - outW(j, i, 1)) ^ 2 + (cityPos(2) - outW(j, i, 2)) ^ 2;
                end
            end
        end
%       找出获胜节点并记录，更新其状态
        [minRow, minCol] = find(dis == min(min(dis)));
        done(minCol(1)) = 1;
    
        neigborFun = zeros(1, outNum);
        for i = 1 : outNum
%           计算环形结构上节点i和获胜节点minIndex的距离
            d = min(abs(i - minCol(1)), outNum - abs(i - minCol(1)));
%           对获胜节点邻域的所有节点进行权值更新
            if d < neigborPercent * outNum
                neigborFun(i) = exp(- d * d / (g * g));
                outW(minRow(1),i,1) = outW(minRow(1),i,1) + beta * neigborFun(i) * (cityPos(1)-outW(minRow(1),i,1));
                outW(minRow(1),i,2) = outW(minRow(1),i,2) + beta * neigborFun(i) * (cityPos(2)-outW(minRow(1),i,2));
            end
        end
        
        % 若获胜节点到城市的距离小于0.01*mindst则权值调整为城市坐标
        if sqrt((cityPos(1) - outW(minRow(1), minCol(1), 1)).^2 + (cityPos(2) - outW(minRow(1),minCol(1), 2)).^2) < 0.01 * mindst
            outW(minRow(1),minCol(1), 1) = cityPos(1);
            outW(minRow(1),minCol(1), 2) = cityPos(2);
%             if cityOver(minCol(1)) == 0
                stopCity = stopCity + 1;
%                 cityOver(minCol(1)) = 1;
%             end
%             resTour(recNum) = minIndex;
%             recNum = recNum + 1;
              preserveNode(minRow(1),countTourCityNum(minRow(1)),1) =  cityPos(1);
              preserveNode(minRow(1),countTourCityNum(minRow(1)),2) =  cityPos(2);
              countTourCityNum(minRow(1)) = countTourCityNum(minRow(1)) + 1;
        end
    end
%   调整减小学习率，是迭代过程逐渐收敛
    alpha = alpha * 0.998;
    g = (1 - alpha) * g;
%   当权值没有改变时，结束训练
    if outW == lastW
        break;
    elseif stopCity >= cityNum
        break;
    end
    lastW = outW;   %重新保存上一次的权值
end
% 计时结束
toc


% 导入最短路径数据
% cityTour = load('eil101.tour.txt');
% cityTour = load('att48.tour.txt');
% cityTour = load('bayg29.tour.txt');
% cityTour = load('eil51.tour.txt');
% cityTour = load('eil76.tour.txt');
cityTour = load('pa561.tour.txt');
% cityTour = load('tsp225.tour.txt');



% 计算最短路径长度(一个环路的）
tourMinDst = 0;
for i = 1 : cityNum - 1
    tourMinDst = tourMinDst + norm(city(cityTour(i + 1)) - city(cityTour(i)), 2);
end
tourMinDst = tourMinDst + norm(city(cityTour(cityNum)) - city(cityTour(1)), 2);
  

% 计算近似解长度（可能为多环路之和）
tourAppDst = zeros(tourNum,1);
for j = 1 : tourNum
    for i = 1 : outNum - 1
        tourAppDst(j) = tourAppDst(j) + norm(outW(j,i + 1) - outW(j,i), 2);
    end
    tourAppDst(j) = tourAppDst(j) + norm(outW(j,outNum) - outW(j,1), 2); 
end
tourSum = sum(tourAppDst);  %若为多环路则计算多环路之和


% saveas(gcf,'test.png');
title(['loss = ',num2str((tourSum - tourMinDst) / tourMinDst * 100),'% dist = ', num2str(tourSum)]);


fprintf('%d %d %.0f %.2f%%',cityNum, tourMinDst, tourSum, (tourSum - tourMinDst) / tourMinDst * 100);
% hold on;
% % 画出最短路径图
% plot(city(:,1), city(:,2), 'or');
% % title(['loss = ',num2str((tourAppDst - tourMinDst) / tourMinDst * 100),'%']);
% title(['the shortest distance = ', num2str(tourMinDst)]);
% hold on;
% plot(city(cityTour,1), city(cityTour,2), '--g');
% plot([city(cityTour(1),1),city(cityTour(cityNum),1)], [city(cityTour(1),2),city(cityTour(cityNum),2)], '--g');
% pause(0.001);
% hold off;






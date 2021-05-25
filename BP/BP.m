% 1 载入数据
clear all;
load data.dat

% 2 样本数据处理

% 打乱样本顺序
randIndex = randperm(size(data,1));
inputData = data(randIndex,:);

% 样本属性通用标准化[0,1]
for i = 3 : 6
    inputData(:,i) = mapminmax(inputData(:,i)',-1,1)';
end

inputData(:,end+1:end+3) = 0;

for i = 1 : size(inputData,1)
    switch inputData(i,2)
        case 0
            inputData(i,end-2:end) = [1,0,0];
        case 1
            inputData(i,end-2:end) = [0,1,0];
        case 2
            inputData(i,end-2:end) = [0,0,1];
        otherwise  
    end
    
end

trainNum = size(inputData,1) * 0.8; %训练样本数
trainData = inputData(1:trainNum,:); %训练样本数据
testData = inputData(trainNum + 1:end,:); %测试样本数据
testNum = size(testData,1); %测试样本数

% 初始化权值[-1,1]
nL1 = 4;  %输入层结点数
nL2 = 3;  %隐含层结点数
nL3 = 3;  %输出层结点数

V = unifrnd(-1,1,nL1,nL2);  %输入层到隐含层的权值
W = unifrnd(-1,1,nL2,nL3);  %隐含层到输出层的权值
doorL2 = unifrnd(-1,1,nL2,1)';  %隐含层阈值
doorL3 = unifrnd(-1,1,nL3,1)';  %输出层阈值

actL2 = [];  %隐藏层激活值
actL3 = [];  %输出层激活值

eL3 = [];  %输出层误差参数
eL2 = []; %隐含层误差参数

Alpha = 0.001;  %学习率1
Beta = 0.001;   %学习率2
% Gamma = 0.8  %过去权值常量

Maxpoch = 1000; %最大迭代次数
finEk = [0]; %记录每次迭代的累计误差
ekStop = 1e-5; %判断误差收敛的预设值

for i = 1 : Maxpoch
    ek = 0;
    for k = 1 : trainNum
        actL2 = extractdata(sigmoid(dlarray(trainData(k,3:6) * V+doorL2)));
        actL3 = extractdata(sigmoid(dlarray(actL2 * W+doorL3)));
        
        ek = ek + sum((trainData(k,end-2:end)-actL3).^2,2)/2;
        eL3 = actL3.*(1-actL3).*(trainData(k,end-2:end)-actL3);
        eL2 = actL2.*(1-actL2).*(eL3*W');
        W = W + actL2'*eL3*Alpha;
        V= V + trainData(k,3:6)'*eL2*Beta;
        doorL2 = doorL2 + eL2*Beta;
        doorL3 = doorL3 + eL3*Alpha;     
    end
    
    ek = ek / trainNum;
    if (abs(ek - finEk(end)) <= ekStop)
        break;
    end
    finEk = [finEk;ek];
end


%绘制损失函数变化曲线
plot(finEk(2:end,:));


%测试
rightNum = 0;
forecast = [];
for i = 1 : testNum
    actL2 = extractdata(sigmoid(dlarray(testData(i,3:6) * V+doorL2)));
    actL3 = extractdata(sigmoid(dlarray(actL2 * W+doorL3)));
    [~,index] = max(actL3,[],2);
    if (index - 1 == testData(i,2))
        rightNum = rightNum + 1;
    end
    forecast = [forecast;index-1];
end

%精确度
accuracy = 100 * rightNum / testNum;
disp([num2str(accuracy) '%']);

%绘制预测与测试差距
N = 1:testNum;
figure; plot(N,forecast,'bo',N,testData(:,2),'r*');

%打印权值和阈值




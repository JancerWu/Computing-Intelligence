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

trainNum = size(inputData,1) * 0.8;
trainData = inputData(1:trainNum,3:end);
testData = inputData(trainNum + 1:end,3:end);
testNum = size(testData,1);

net = newrb(trainData(:,1:4)',trainData(:,end - 2: end)',0.0,0.6);
res = net(testData(:,1:4)')';

rightNum = 0;
for i = 1 : size(testData,1)
    [~,index] = max(res(i,:),[],2);
    if (index == (find(testData(i,:)==1) - 4))
        rightNum = rightNum + 1;
    end
end

%精确度
accuracy = 100 * rightNum / testNum;
disp([num2str(accuracy) '%']);

function [mindst] = findmin_distance(datanew)
distance = zeros((size(datanew)).^2);
% 找到二维点中的最小距离
sizedata = size(datanew);
k = 1;
mindst = +inf;
for i=1:sizedata
    for j=i:sizedata
        distance(k) = sqrt((datanew(i,1)-datanew(j,1)).^2+ (datanew(i,2)-datanew(j,2)).^2);
        k = k +1;
    end
end
for i=1:size(distance)
    if(distance(i)~=0 && distance(i) < mindst)
       mindst =  distance(i);
    end
end
end


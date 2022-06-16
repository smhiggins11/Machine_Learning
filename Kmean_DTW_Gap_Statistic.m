
Log_Wk = log(Within_Cluster_Total_Sum);

%% 

Matrix_Sum = zeros(101,1);

for ii = 1:size(Matrix,2)
    Matrix_Sum = Matrix_Sum+Matrix(:,ii);
end

Matrix_Ave = Matrix_Sum/size(Matrix,2);

for ii = 1:size(Matrix,2)
    Matrix_Deviation(:,ii) = (Matrix(:,ii)-Matrix_Ave).^2;
end

Matrix_Deviation_Sum = zeros(101,1);
for ii = 1:size(Matrix,2)
    Matrix_Deviation_Sum = Matrix_Deviation_Sum+Matrix_Deviation(:,ii);
end

Matrix_Deviation_Ave = Matrix_Deviation_Sum/size(Matrix,2);

Matrix_SD_plus = sqrt(Matrix_Deviation_Ave);
Matrix_SD_plus_Ave = mean(Matrix_SD_plus);

Matrix_SD_min = Matrix_SD_plus*-1;
Matrix_SD_min_Ave = mean(Matrix_SD_min);

% for ii = 1:size(Matrix,2)
%     for b = 1:length(Matrix)
%         Matrix_Max(b,ii) = Matrix(b,ii)+Matrix_SD_plus(b,1);
%     end
% end

for ii = 1:size(Matrix,2)
    for b = 1:length(Matrix)
        Matrix_Max(b,ii) = Matrix(b,ii)+Matrix_SD_plus_Ave;
    end
end

% for ii = 1:size(Matrix,2)
%     for b = 1:length(Matrix)
%         Matrix_Min(b,ii) = Matrix(b,ii)+Matrix_SD_min(b,1);
%     end
% end

for ii = 1:size(Matrix,2)
    for b = 1:length(Matrix)
        Matrix_Min(b,ii) = Matrix(b,ii)+Matrix_SD_min_Ave;
    end
end

plot(Matrix_Max)

Elog_Wk = zeros(1,100);
% Data_Max = zeros(length(Matrix),1);
% for ii = 1:length(Matrix)
%     Data_Max(ii,1) = max(round(Matrix(ii,:)));
% end
% 
% Data_Min = zeros(length(Matrix),1);
% for ii = 1:length(Matrix)
%     Data_Min(ii,1) = min(round(Matrix(ii,:)));
% end

N = 1;
while N ~=101    
% Generated_Data = zeros(length(Matrix),size(Matrix,2));
% for ii = 1:size(Matrix,2)
%     for b = 1:length(Matrix)
%         Generated_Data(b,ii) = Matrix_Max(b,ii) + (Matrix_Min(b,ii)-Matrix_Max(b,ii)) .* rand(1,1);
%     end
% end

Generated_Data = zeros(length(Matrix),size(Matrix,2));
for ii = 1:size(Matrix,2)
    Generated_Data(:,ii) = Matrix(:,ii)+(Matrix_SD_plus_Ave + (Matrix_SD_min_Ave-Matrix_SD_plus_Ave) .* rand(1,1));
end

Cluster_Centroid = zeros(length(Matrix),k);
for ii = 1:k
    Cluster_Centroid(:,ii) = Generated_Data(:,Cluster_Rand(ii));
end
plot(Generated_Data)
%% K-means on generated data
Cluster_Converge_Ave = 1;%I chose one for this variable as an arbitrary number just to create the variable before the while loop.
n = 0; %This variable will tell me how many time k-means runs before it stops.
while Cluster_Converge_Ave ~= 0 
Cluster_DTW = zeros(1,size(Generated_Data,2));
for ii = 1:k
    Cluster_DTW(ii,:) = zeros(1,size(Generated_Data,2));
end

%determines the distance (using dtw) between the initial centroids and all the vectors in the matrix 
for ii = 1:k
    for a=1:size(Generated_Data,2)
    Cluster_DTW(ii,a) = dtw(Cluster_Centroid(:,ii), Generated_Data(:,a));
    end
end

%in order for me to group all the Cluster_DTW variables together, I need to round the data.
for ii = 1:k
   Cluster_DTW(ii,:) = round(Cluster_DTW(ii,:));
end

%Here I find which centroid each vector is closest to by finding the
%minimum DTW distance from all three centroids
Cluster_min = double(1:size(Generated_Data,2));
for ii = 1:size(Generated_Data,2)
    Cluster_min(ii) = min(Cluster_DTW(:,ii));
end

%Now I use the minimum distance found above and will match that distance
%measure to the appropriate vector
Cluster_Num = zeros(k,size(Generated_Data,2));
Cluster_length = zeros(1,k);
Cluster_Num_1 = cell(1,k);
for ii = 1:k
    Cluster_Num(ii,:) = zeros(1,size(Generated_Data,2));
    Cluster_length(ii) = length(find(Cluster_DTW(ii,:)==Cluster_min));
    Cluster_Num(ii,1:Cluster_length(1,ii)) = find(Cluster_DTW(ii,:)==Cluster_min);
    Cluster_Num_1{ii} = {Cluster_Num(ii,:)};
end

%removes zeros from cell arrays
for ii = 1:k
    Cluster_Num_1{1,ii}{1,1}(Cluster_Num_1{1,ii}{1,1} == 0) = [];
end

%Here I will take those find the vectors for each cluster within the Matrix
Cluster_Vectors_Cell = cell(1,k);
for ii = 1:k
    Cluster_Vectors_Doub = Generated_Data(:,Cluster_Num_1{1,ii}{1,:});
    Cluster_Vectors_Cell{ii} = {Cluster_Vectors_Doub};
    Cluster_Vectors_Doub = [];
end

%Sum all the vectors within each Cluster
Cluster_Sum = cell(1,k);
for ii = 1:k
    Cluster_Sum{ii} = 0;
    for a = 1:length(Cluster_Num_1{1,ii}{1,1})
        Cluster_Sum{ii} = Cluster_Sum{ii} + Cluster_Vectors_Cell{1,ii}{1,1}(:,a);
    end
end

%Average each cluster
Cluster_Ave = cell(1,k);
for ii = 1:k
    Cluster_Ave{ii} = Cluster_Sum{ii}/length(Cluster_Num_1{1,ii}{1,1});
end

%Find how many participants are in each cluster
for ii = 1:k
    Cluster_length(ii) = length(Cluster_Num_1{1,ii}{1,1}(1,:));
end

%This will compare the subject numbers from the previous run of k-means to
%the current run of k-means to determine if the clusters changed. If the
%subjects within each cluster do not change, k-means will stop.
%This for loop will only start after the first complete cycle of the while
%loop.
if n > 0
    Cluster_Num_DTW = cell(k,size(Generated_Data,2));
    Cluster_Num_Difference = cell(1,k);
    for ii = 1:k
        for a = 1:k
        Cluster_Num_DTW{ii,a} = dtw(Cluster_Num_1{1,ii}{1,1},Cluster_Num_Previous{1,a}{1,1});
        if Cluster_Num_DTW{ii,a} == 0
            Cluster_Num_Difference{ii} = Cluster_Num_DTW{ii,a};
        end
        end
    end 
Cluster_Num_DTW_Double = cell2mat(Cluster_Num_DTW);
Cluster_Converge = zeros(1,k);
for ii = 1:k
    Cluster_Converge(ii) = Cluster_Num_DTW_Double(ii,ii);
end
Cluster_Converge_Ave = mean(Cluster_Converge);
end

%change variable with subject numbers to compare to the next run of k-means
Cluster_Num_Previous = Cluster_Num_1; 

%I now change Cluster average variables to cluster centroid variables in
%order to replace the intial centroids at the beginning of the while loop
for ii = 1:k
    Cluster_Centroid(:,ii) = Cluster_Ave{1,ii}(:,1);
end
n = n+1;
end

idx = zeros(size(Generated_Data,2),1);
for ii = 1:size(Cluster_Num_1,2)
    for a = 1:length(Cluster_Num_1{1,ii}{1,1})
        idx(Cluster_Num_1{1,ii}{1,1}(1,a),1) = ii;
    end
end
for ii = 1:k
    for b = 1:size(Cluster_Vectors_Cell{1,ii}{1,1},2)
        Gen_Within_Cluster(ii).Distance(1,b) = dtw(Cluster_Vectors_Cell{1,ii}{1,1}(:,b), Cluster_Centroid(:,ii));
    end
end

for ii = 1:k
    Gen_Within_Cluster_Sum(1,ii) = sum(Gen_Within_Cluster(ii).Distance);
end

Gen_Within_Cluster_Total_Sum = sum(Gen_Within_Cluster_Sum);

Elog_Wk(1,N) = log(Gen_Within_Cluster_Total_Sum);
N = N+1;
end

Elog_Wk_Ave = mean(Elog_Wk);

Gap_Stat = Elog_Wk_Ave - Log_Wk;
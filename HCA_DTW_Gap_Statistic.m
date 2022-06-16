%August 26 2021
%Gap statistic calculation for HCA DTW Ward
clear;
clc;

filename =('Pelvis_Flex_Ext_Kinematics_Data_Final.xlsx');%load kinematic file
Matrix = xlsread(filename);

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

for ii = 1:size(Matrix,2)
    for b = 1:length(Matrix)
        Matrix_Max(b,ii) = Matrix(b,ii)+Matrix_SD_plus_Ave;
    end
end

for ii = 1:size(Matrix,2)
    for b = 1:length(Matrix)
        Matrix_Min(b,ii) = Matrix(b,ii)+Matrix_SD_min_Ave;
    end
end

plot(Matrix_Max)

Elog_Wk = zeros(1,100);

N = 1;
while N ~=101  
    Generated_Data = zeros(length(Matrix),size(Matrix,2));
for ii = 1:size(Matrix,2)
    Generated_Data(:,ii) = Matrix(:,ii)+(Matrix_SD_plus_Ave + (Matrix_SD_min_Ave-Matrix_SD_plus_Ave) .* rand(1,1));
end

Generated_Data_Sum = Generated_Data;
Generated_Data1 = Generated_Data;

w = 1;

Clusters = cell(1,size(Generated_Data,2));
for k = 1:size(Generated_Data1,2)
    Clusters{k} = k;
end

plot(Generated_Data)

%% HCA calculation
Number_of_Clusters = 0;%This variable is used to control how many clusters are created

%Preallocation of data that will be used to create the dendrogram
Dendrogram_Total = zeros(size(Generated_Data,2)-1,3);

%while loop is used to run throught the clustering algorithm until a number
%of clusters specific has been reached. 
while Number_of_Clusters ~=10 %change this variable to change the number of clusters produced
    
%The while loop and for loop below will calculate the dtw between the first
%vector and all the other vectors in the matrix. The while loop will then
%switch to compare the distance between the 2nd vector and the other vectors
%in the matrix and this will be done for each vector in the matrix.
SSE_Dist = cell(1,size(Generated_Data,2));%preallocation of variable
n = 1;%start n at 1
while n < size(Generated_Data,2)+1
    for k = 1:size(Generated_Data,2)
        if Clusters{1,n} == 0
            SSE_Dist{1,n}(ii,k) = 0;
        elseif Clusters{1,n} ~= 0
            for ii = 1:length(Clusters{1,n})
                SSE_Dist{1,n}(ii,k) = dtw(Generated_Data(:,n),Generated_Data(:,k));
            end
        end
    end
n = n+1;
end

%calculates the square root of 2*(number of subjects in cluster 1)*(number of subjects in cluster 2)
%divided by (number of subjects in cluster 1)+(number of subjects in cluster 2)
SSE_sqrt = cell(1,size(Generated_Data,2));
n = 1;
while n < size(Generated_Data,2)+1
    for k = 1:size(Generated_Data,2)
        SSE_sqrt{1,n}(1,k) = sqrt(2*(size(Clusters{n},2))*(size(Clusters{k},2))/((size(Clusters{n},2))+(size(Clusters{k},2))));
    end
    n = n+1;
end

%creates a matrix through multiplying the SSE_sqrt variable by the SSE_Dist
%variable
SSE = zeros(size(Generated_Data,2),size(Generated_Data,2));
n = 1;
while n < size(Generated_Data,2)+1
    for k = 1:size(Generated_Data,2)
      SSE(n,k) = SSE_sqrt{1,n}(1,k)*SSE_Dist{1,n}(1,k);  
    end
    n = n+1;
end

%rearranges the SSE variable to place all the zeros in the first column
SSE_Rearranged = zeros(size(Generated_Data,2),size(Generated_Data,2));
for k = 1:size(Generated_Data,2)
    SSE_Rearranged(k,:) = sort(SSE(k,:));
end

%deletes the first column containing all zeros
SSE_Rearranged(:,1) = [];

%Now the first column contains all the minimum SSE values for each
%time-seris
SSE_Min_Column = SSE_Rearranged(:,1);

%as the main while loop continues to run through this part of the
%SSE_Min_Column will start to accumulate zero values. This section will
%change those zeros to NaN.
for ii = 1:length(SSE_Min_Column)
    if SSE_Min_Column(ii,1) == 0
        SSE_Min_Column(ii,1) = NaN;
    end
end

%finds the minimum value in the SSE_Min_Column variable. 
SSE_Min = min(SSE_Min_Column);

%This will find the participant column numbers that match the SSE_Min value
SSE_Min_Find = find(SSE_Min_Column == SSE_Min);

%The Clusters variable is used to organize the participants into clusters
%This for loop will place the participants numbers from the SSE_Min_Find
%variable and places those numbers at the end of the Clusters variable. I
%need two for loops to place those numbers one at a time. 
for ii = 1:length(Clusters{1,SSE_Min_Find(1,1)}) 
    Clusters{:,size(Generated_Data,2)+1}(1,ii) = Clusters{1,SSE_Min_Find(1,1)}(1,ii);      
end

for ii = 1:length(Clusters{1,SSE_Min_Find(2,1)}) 
    Clusters{:,size(Generated_Data,2)+1}(1,end+1) = Clusters{1,SSE_Min_Find(2,1)}(1,ii);      
end

%sums the time-seris for the participants from othe SSE_Min_Find variable
Cluster_Sum = zeros(101,1);
for ii = 1:length(SSE_Min_Find)
    Cluster_Sum = Cluster_Sum + Generated_Data_Sum(:,SSE_Min_Find(ii,1));
end

%This adds the sum of the vectors at the end of the Matrix_Sum variable
Generated_Data_Sum(:,size(Generated_Data_Sum,2)+1) = Cluster_Sum;

%I do not want to reuse the time-series again after each time the main
%while loop runs. Therefore, for the participants combined to create a new
%cluster, I will change the original columns containing their data to zeros
for ii = 1:length(SSE_Min_Find)
   Generated_Data_Sum(:,SSE_Min_Find(ii,1)) = 0;%Need to change the value to something else that works.
end

%This variable finds the time-series that match the columns from the
%SSE_Min_Find variable
Hierarchical_Vector_find = zeros(size(Generated_Data,1),length(Clusters{end}));
for ii = 1:length(Clusters{end})
    Hierarchical_Vector_find(:,ii) = Generated_Data1(:,Clusters{end}(1,ii));
end

%this sums the time-seris from Hierarchical_Vector_Find variable
Hierarchical_Vector_Sum = zeros(101,1);
for ii = 1:size(Hierarchical_Vector_find,2)
    Hierarchical_Vector_Sum = Hierarchical_Vector_Sum(:,1)+Hierarchical_Vector_find(:,ii);
end    

%this averages the Hierarchical_Vector_Sum variable
Hierarchical_Vector_Ave = Hierarchical_Vector_Sum/size(Hierarchical_Vector_find,2);

%Change columns from SSE_Min_Find to a new value. 
%This is used to eliminate those column but also keeping the column numbers
for ii = 1:length(SSE_Min_Find)
    Generated_Data(:,SSE_Min_Find(ii,1)) = -11111111111111111111111111;%Need to change the value to something else that works.
end

%This adds the new averaged cluster to the end of the Matrix data
Generated_Data(:,size(Generated_Data,2)+1) = Hierarchical_Vector_Ave;

%Now I need to place the two clusters being combined and the distance
%between those two clusters in the Dendrogram variable. This variable will
%be updata after each run of the main while loop
Dendrogram = [];
Dendrogram(1,:) = SSE_Min_Find;
Dendrogram(1,3) = SSE_Min;

Dendrogram_Total(w,:) = Dendrogram;

%adds a zeros to the columns in the Clusters variable that have already
%been used.
for ii = 1:length(SSE_Min_Find)
    Clusters{SSE_Min_Find(ii,1)} = 0;
end

%this finds which column contain clusters (i.e. columns that don't equal zero)
find_clusters = cell(1,1);
for ii = 1:length(Clusters)
    find_clusters{ii} = find(Clusters{ii}(1,1) > 0);
end

%This finds the empty cell and eliminates them
emptyCells = cellfun('isempty', find_clusters);
find_clusters(emptyCells) = [];

%This varible contains how man clusters there
Number_of_Clusters = length(find_clusters);

w = w+1;%w will increase by 1 after each run of the main while loop
end

%All cells that contain a zero will be changed to empty
for ii = 1:length(Clusters)
    if Clusters{ii} == 0
        Clusters{ii} = [];
    end
end

%This will find all empty cell and delete them leaving only the cells that
%contain the participants in each cluster
emptyCells = cellfun('isempty', Clusters);
Clusters(emptyCells) = [];

%Finds all the time-seris within each cluster
Cluster_Vectors = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_Vectors{ii} = Generated_Data1(:,Clusters{ii});
end

%Sums up all the time-series within each cluster
Clusters_Sum = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Clusters_Sum{ii} = zeros(length(Generated_Data1),1);
end

for ii = 1:Number_of_Clusters
    for a = 1:size(Cluster_Vectors{ii},2)
        Clusters_Sum{ii} = Clusters_Sum{ii} + Cluster_Vectors{ii}(:,a);
    end
end

%averages each cluster
Cluster_Ave = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_Ave{ii} = Clusters_Sum{ii}/size(Cluster_Vectors{ii},2);
end

%% Total Within Sum of Squares
Within_Cluster = struct('Distance',cell(1,Number_of_Clusters));%preallocate variable
for ii = 1:Number_of_Clusters
    for a = 1:size(Cluster_Vectors{1,ii},2)
        Within_Cluster(ii).Distance(a,1) = dtw(Cluster_Vectors{1,ii}(:,a), Cluster_Ave{:,ii}(:,1));
    end
end

Within_Cluster_Sum = zeros(1,Number_of_Clusters);%preallocate variable
for ii = 1:Number_of_Clusters
    Within_Cluster_Sum(ii) = sum(Within_Cluster(ii).Distance);
end

Within_Cluster_Total_Sum = sum(Within_Cluster_Sum);
Elog_Wk(1,N) = log(Within_Cluster_Total_Sum);
N = N+1;
end
Elog_Wk_Ave = mean(Elog_Wk);

%Gap_Stat = Elog_Wk_Ave - Log_Wk;

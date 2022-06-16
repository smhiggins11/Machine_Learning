clear;
clc;
%% Upload all files and create each variable
filename = ('Lumbar_Flex_Ext_Kinematics_Data_Final.xlsx');%load file containing kinematics
filename1 = ('Norm_Running_Age_Sex_Height_Weight.xlsx');%load file containing anthropometrics

%I created three identical matrix variables that will be used differently
%throughout the code
Matrix = xlsread(filename);%Average of the new clusters will be added to this variable
Matrix1 = xlsread(filename);%This variable will not change
Matrix_Sum = xlsread(filename);%Sum of the new clusters will be added to this variable

%Here are the different anthropometric data I will use after cluster has
%been completed
[num,text] = xlsread(filename1);
Subjects = text(2:length(text),1);%subject ID
Gender = text(2:length(text),2);
Age = num(:,1);
Height = num(:,2);
Weight = num(:,3);

w = 1;%variable tells me how many times the while loop runs

%Creates a 1xN cell containing subject cell numbers
%this will be used to organize subject numbers in clusters
Clusters = cell(1,size(Matrix,2));
for k = 1:size(Matrix1,2)
    Clusters{k} = k;
end

%variable starts at 0 to create the variable before starting the while loop
Number_of_Clusters = 0;%This variable is used to control how many clusters are created

%Preallocation of data that will be used to create the dendrogram
Dendrogram_Total = zeros(size(Matrix,2)-1,3);

%while loop is used to run throught the clustering algorithm until a number
%of clusters specific has been reached. 
while Number_of_Clusters ~=4 %change this variable to change the number of clusters produced
    
%The while loop and for loop below will calculate the dtw between the first
%vector and all the other vectors in the matrix. The while loop will then
%switch to compare the distance between the 2nd vector and the other vectors
%in the matrix and this will be done for each vector in the matrix.
SSE_Dist = cell(1,size(Matrix,2));%preallocation of variable
n = 1;%start n at 1
while n < size(Matrix,2)+1
    for k = 1:size(Matrix,2)
        if Clusters{1,n} == 0
            SSE_Dist{1,n}(ii,k) = 0;
        elseif Clusters{1,n} ~= 0
            for ii = 1:length(Clusters{1,n})
                SSE_Dist{1,n}(ii,k) = dtw(Matrix(:,n),Matrix(:,k));
            end
        end
    end
n = n+1;
end

%calculates the square root of 2*(number of subjects in cluster 1)*(number of subjects in cluster 2)
%divided by (number of subjects in cluster 1)+(number of subjects in cluster 2)
SSE_sqrt = cell(1,size(Matrix,2));
n = 1;
while n < size(Matrix,2)+1
    for k = 1:size(Matrix,2)
        SSE_sqrt{1,n}(1,k) = sqrt(2*(size(Clusters{n},2))*(size(Clusters{k},2))/((size(Clusters{n},2))+(size(Clusters{k},2))));
    end
    n = n+1;
end

%creates a matrix through multiplying the SSE_sqrt variable by the SSE_Dist
%variable
SSE = zeros(size(Matrix,2),size(Matrix,2));
n = 1;
while n < size(Matrix,2)+1
    for k = 1:size(Matrix,2)
      SSE(n,k) = SSE_sqrt{1,n}(1,k)*SSE_Dist{1,n}(1,k);  
    end
    n = n+1;
end

%rearranges the SSE variable to place all the zeros in the first column
SSE_Rearranged = zeros(size(Matrix,2),size(Matrix,2));
for k = 1:size(Matrix,2)
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
    Clusters{:,size(Matrix,2)+1}(1,ii) = Clusters{1,SSE_Min_Find(1,1)}(1,ii);      
end

for ii = 1:length(Clusters{1,SSE_Min_Find(2,1)}) 
    Clusters{:,size(Matrix,2)+1}(1,end+1) = Clusters{1,SSE_Min_Find(2,1)}(1,ii);      
end

%sums the time-seris for the participants from othe SSE_Min_Find variable
Cluster_Sum = zeros(101,1);
for ii = 1:length(SSE_Min_Find)
    Cluster_Sum = Cluster_Sum + Matrix_Sum(:,SSE_Min_Find(ii,1));
end

%This adds the sum of the vectors at the end of the Matrix_Sum variable
Matrix_Sum(:,size(Matrix_Sum,2)+1) = Cluster_Sum;

%I do not want to reuse the time-series again after each time the main
%while loop runs. Therefore, for the participants combined to create a new
%cluster, I will change the original columns containing their data to zeros
for ii = 1:length(SSE_Min_Find)
    Matrix_Sum(:,SSE_Min_Find(ii,1)) = 0;%Need to change the value to something else that works.
end

%This variable finds the time-series that match the columns from the
%SSE_Min_Find variable
Hierarchical_Vector_find = zeros(size(Matrix,1),length(Clusters{end}));
for ii = 1:length(Clusters{end})
    Hierarchical_Vector_find(:,ii) = Matrix1(:,Clusters{end}(1,ii));
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
    Matrix(:,SSE_Min_Find(ii,1)) = -11111111111111111111111111;%Need to change the value to something else that works.
end

%This adds the new averaged cluster to the end of the Matrix data
Matrix(:,size(Matrix,2)+1) = Hierarchical_Vector_Ave;

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
    Cluster_Vectors{ii} = Matrix1(:,Clusters{ii});
end

%Sums up all the time-series within each cluster
Clusters_Sum = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Clusters_Sum{ii} = zeros(length(Matrix1),1);
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

%Creates a vector containing which cluster number each participant belongs
%too. This will be used to calculate the dunn index.
idx = zeros(1,size(Matrix1,2));
for ii = 1:size(Clusters,2)
    for a = 1:length(Clusters{1,ii})
        idx(1,Clusters{1,ii}(1,a)) = ii;
    end
end

%This variable caclulates the dtw distance between each vector and all the
%other vectors in the matrix. This will be used to cacluate the dunn index.
distM = zeros(size(Matrix1,2),size(Matrix1,2));
n = 1;
while n < size(Matrix1,2)+1
    for k = 1:size(Matrix1,2)
        distM(n,k) = dtw(Matrix1(:,n),Matrix1(:,k)); 
    end
    n = n+1;
end

%Calculates dunn index; dunn(number of clusters,distM,indx)
Dunn = dunns(Number_of_Clusters,distM,idx);

%Creates the dendrogram only when the Number_of_Cluster equal 1
if Number_of_Clusters == 1
    leaforder = optimalleaforder(Dendrogram_Total, distM);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    H = dendrogram(Dendrogram_Total,0,'Reorder',leaforder,'ColorThreshold',0.5*max(Dendrogram_Total(:,3)),'Reorder',leaforder);
    set(H,'LineWidth',2.0)
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName','Times New Roman','fontsize',15,'FontWeight','bold')
    title('Dendrogram For Lumbar Kinematics', 'FontName','Times New Roman','FontSize', 28,'FontWeight','bold')
    ylabel('Ward Linkage Distance','FontName','Times New Roman','FontSize',22,'FontWeight','bold')
    xlabel('Individual Subjects','FontName','Times New Roman','FontSize',22,'FontWeight','bold')
%     for ii = 1:47
%         H(ii).Color = 'b';
%     end
set(gcf,'color','w');%changes graph background from grey to white
end
   
%pairs the participant ID's with the participant column numbers
Cluster_Subjects = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_Subjects{1,ii} = Subjects(Clusters{1,ii});
end

%finds the ages of each participant within each cluster
Cluster_Age = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_Age{1,ii} = Age(Clusters{1,ii});
end

%finds the average age within each cluster
Cluster_Age_Ave = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_Age_Ave{1,ii} = mean(Cluster_Age{1,ii});
end

%finds the height of each participant within each cluster
Cluster_Height = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_Height{1,ii} = Height(Clusters{1,ii});
end

%averages the height within each cluster
Cluster_Height_Ave = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_Height_Ave{1,ii} = mean(Cluster_Height{1,ii});
end

%finds the weight of each participant within each cluster
Cluster_Weight = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_Weight{1,ii} = Weight(Clusters{1,ii});
end

%finds the average weight within each cluster
Cluster_Weight_Ave = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_Weight_Ave{1,ii} = mean(Cluster_Weight{1,ii});
end

%Finds the number of participants within each cluster
Cluster_length = zeros(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_length(ii) = length(Clusters{ii});
end

Cluster_YA = cell(1,Number_of_Clusters);
Cluster_MA = cell(1,Number_of_Clusters);
for ii = 1:Number_of_Clusters
    Cluster_YA{1,ii} = find(Cluster_Age{1,ii}<41);
    Cluster_MA{1,ii} = find(Cluster_Age{1,ii}>40);
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
%% 
%creates variable used to display cluster number and number of participants
%in each cluster when I plot the final clustering result.
Cluster_Total = zeros(1,Number_of_Clusters*2);
a = 0;
for ii = 1:Number_of_Clusters
    Cluster_Total(ii+a) = ii;
    a = a+1;
    Cluster_Total(ii+a) = (Cluster_length(ii));
end

%figure
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
hold on
for ii = 1:Number_of_Clusters
    plot(Cluster_Ave{:,ii},'LineWidth',4.0)
end
hold on
set(gca,'XTick',[0 10 20 30 40 50 60 70 80 90 100],'FontName','Times New Roman','FontSize',18,'FontWeight','bold')
xlim([0 100])
set(gca,'YTick',[-50 -40 -30 -20 -10 0 10 20],'FontName','Times New Roman','FontSize',18,'FontWeight','bold')
ylim([-50 20])
xline(34,'--');
xline(54,'--');
xline(73,'--');
title('HCA Clustering For Pelvis Kinematics','FontName','Times New Roman','FontSize', 28,'FontWeight','bold')
set(get(gca,'title'),'Position',[50 22 1.00011])
xlabel('Trunk Flexion and Extension Cycle (%)','FontName','Times New Roman','FontSize',22,'FontWeight','bold')
set(get(gca,'xlabel'),'Position',[50 -55 1.00011])
ylabel('Trunk Flexion/Extension Angle (°)','FontName','Times New Roman','FontSize',22,'FontWeight','bold')
hold on
annotation('textbox',[.11 .070 .2 .01],'String','Neutral','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.36 .070 .2 .01],'String','Max Flexion','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.53 .070 .2 .01],'String','Neutral','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.66 .070 .2 .01],'String','Hyperextension','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.885 .070 .2 .01],'String','Neutral','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
legend(sprintfc('Cluster %d (n = %d)',Cluster_Total),'Location','east','FontName','Times New Roman','FontSize', 18,'FontWeight','bold')
set(gcf,'color','w');%changes graph background from grey to white
legend boxoff
clear;
clc;

filename =('Lumbar_Flex_Ext_Kinematics_Data_Final.xlsx');%load kinematic file
filename1 = ('Norm_Running_Age_Sex_Height_Weight.xlsx');%load file containing anthropometrics
Matrix = xlsread(filename);%create time-series data as a matrix.
%Data is nXm matrix (n = time series; m = participant)

%Here are the different anthropometric data I will use after cluster has
%been completed
[num,text] = xlsread(filename1);
Subjects = text(2:length(text),1);%subject ID
Gender = text(2:length(text),2);
Age = num(:,1);
Height = num(:,2);
Weight = num(:,3);

%k = the number of clusters
%change k variable to change the number of clusters calculated
k = 4;

%Generates a random number within the size of the matrix for each cluster
%These random numbers will be used as the initial centroids for each
%cluster. The random number generator is set up to select nonrepeating
%numbers
Cluster_Rand = randperm(size(Matrix,2));
Cluster_Rand = Cluster_Rand(1:k);

%if I want to manually enter cluster centroids do it here. Keep code
%commented out when not using. I only use this when I want to view cluster
%results using centroids I already calculated.
%for pelvis data I used [34 10]
%for lumbar data I used [70 60 3 68]
Cluster_Rand = [70 60 3 68];

%Here I find the time-series that matches the random numbers
%These time-series will be used as the initial centroids to start off
%k-means analysis
Cluster_Centroid = zeros(length(Matrix),k);
for ii = 1:k
    Cluster_Centroid(:,ii) = Matrix(:,Cluster_Rand(ii));
end

%% K-means Algorithm
%This while loop will run k-means clustering until no time-series is reassigned to a different cluster.
%when the "Cluster_Converge_Ave" variable is 0 (no difference in clusters between previous and current interations of k-means)
%k-means analysis will stop. 
Cluster_Converge_Ave = 1;%I chose one for this variable as an arbitrary number just to create the variable before the while loop starts.
n = 0; %This variable will tell me how many time k-means runs before it stops.
while Cluster_Converge_Ave ~= 0 

%Creates kxn matrix of zeros   
Cluster_DTW = zeros(k,size(Matrix,2));%preallocate variable

%determines the distance (using dtw) between the initial centroids and all the vectors in the matrix 
for ii = 1:k
    for a=1:size(Matrix,2)
    Cluster_DTW(ii,a) = dtw(Cluster_Centroid(:,ii),Matrix(:,a));
    end
end

%in order for me to group all the Cluster_DTW variables together, I need to round the data.
for ii = 1:k
   Cluster_DTW(ii,:) = round(Cluster_DTW(ii,:));
end

%Here I find which centroid each vector is closest to by finding the
%minimum DTW distance from all three centroids
Cluster_min = double(1:size(Matrix,2));%preallocate variable
for ii = 1:size(Matrix,2)
    Cluster_min(ii) = min(Cluster_DTW(:,ii));
end

%Now I use the minimum distance found above and will match that distance
%measure to the appropriate vector
Cluster_Num = zeros(k,size(Matrix,2));%preallocate variable
Cluster_length = zeros(1,k);%preallocate variable
Cluster_Num_Cell = cell(1,k);%preallocate variable
for ii = 1:k
    Cluster_Num(ii,:) = zeros(1,size(Matrix,2));
    Cluster_length(ii) = length(find(Cluster_DTW(ii,:)==Cluster_min));
    Cluster_Num(ii,1:Cluster_length(1,ii)) = find(Cluster_DTW(ii,:)==Cluster_min);
    Cluster_Num_Cell{1,ii} = {Cluster_Num(ii,:)};
end

%removes zeros from cell arrays
for ii = 1:k
    Cluster_Num_Cell{1,ii}{1,1}(Cluster_Num_Cell{1,ii}{1,1} == 0) = [];
end

%Here I will find the vectors for each cluster within the Matrix
Cluster_Vectors_Cell = cell(1,k);%preallocate variable
for ii = 1:k
    Cluster_Vectors_Doub = Matrix(:,Cluster_Num_Cell{1,ii}{1,:});
    Cluster_Vectors_Cell{ii} = {Cluster_Vectors_Doub};
    Cluster_Vectors_Doub = [];
end

%Sum all the vectors within each Cluster 
Cluster_Sum = cell(1,k);%preallocate variable
for ii = 1:k
    Cluster_Sum{ii} = 0;
    for a = 1:length(Cluster_Num_Cell{1,ii}{1,1})
        Cluster_Sum{ii} = Cluster_Sum{ii} + Cluster_Vectors_Cell{1,ii}{1,1}(:,a);
    end
end

%Average each cluster
Cluster_Ave = cell(1,k);%preallocate variable
for ii = 1:k
    Cluster_Ave{ii} = Cluster_Sum{ii}/length(Cluster_Num_Cell{1,ii}{1,1});
end

%Find how many participants are in each cluster
for ii = 1:k
    Cluster_length(ii) = length(Cluster_Num_Cell{1,ii}{1,1}(1,:));
end

%This will compare the subject numbers from the previous run of k-means to
%the current run of k-means to determine if the clusters changed. If the
%subjects within each cluster do not change, k-means will stop.
%This for loop will only start after the first complete cycle of the while
%loop.
if n > 0
    Cluster_Num_DTW = cell(k,k);%preallocate variable
    Cluster_Num_Difference = cell(1,k);%preallocate variable
    for ii = 1:k
        for a = 1:k
        Cluster_Num_DTW{ii,a} = dtw(Cluster_Num_Cell{1,ii}{1,1},Cluster_Num_Previous{1,a}{1,1});
        if Cluster_Num_DTW{ii,a} == 0
            Cluster_Num_Difference{ii} = Cluster_Num_DTW{ii,a};
        end
        end
    end 
Cluster_Num_DTW_Double = cell2mat(Cluster_Num_DTW);
Cluster_Converge = zeros(1,k);%preallocate variable
for ii = 1:k
    Cluster_Converge(ii) = Cluster_Num_DTW_Double(ii,ii);
end
Cluster_Converge_Ave = mean(Cluster_Converge);
end

%change variable with subject numbers to compare to the next run of k-means
Cluster_Num_Previous = Cluster_Num_Cell; 

%I now change Cluster average variables to cluster centroid variables in
%order to replace the intial centroids at the beginning of the while loop
for ii = 1:k
    Cluster_Centroid(:,ii) = Cluster_Ave{1,ii}(:,1);
end
n = n+1;
if n == 1000 %if statment will tell me when the while loop has completed 1000 cycle. At this point, the code may not converge and I may need to restart the code.
    disp('Restart')
end
end

%% 
Cluster_Subjects = cell(1,k);
for ii = 1:k
    Cluster_Subjects{1,ii} = Subjects(Cluster_Num_Cell{1,ii}{1,1});
end

%finds the ages of each participant within each cluster
Cluster_Age = cell(1,k);
for ii = 1:k
    Cluster_Age{1,ii} = Age(Cluster_Num_Cell{1,ii}{1,1});
end

%finds the average age within each cluster
Cluster_Age_Ave = cell(1,k);
for ii = 1:k
    Cluster_Age_Ave{1,ii} = mean(Cluster_Age{1,ii});
end

%finds the height of each participant within each cluster
Cluster_Height = cell(1,k);
for ii = 1:k
    Cluster_Height{1,ii} = Height(Cluster_Num_Cell{1,ii}{1,1});
end

%averages the height within each cluster
Cluster_Height_Ave = cell(1,k);
for ii = 1:k
    Cluster_Height_Ave{1,ii} = mean(Cluster_Height{1,ii});
end

%finds the weight of each participant within each cluster
Cluster_Weight = cell(1,k);
for ii = 1:k
    Cluster_Weight{1,ii} = Weight(Cluster_Num_Cell{1,ii}{1,1});
end

%finds the average weight within each cluster
Cluster_Weight_Ave = cell(1,k);
for ii = 1:k
    Cluster_Weight_Ave{1,ii} = mean(Cluster_Weight{1,ii});
end

%Finds the number of participants within each cluster
Cluster_length = zeros(1,k);
for ii = 1:k
    Cluster_length(ii) = length(Cluster_Num_Cell{1,ii}{1,1});
end

Cluster_YA = cell(1,k);
Cluster_MA = cell(1,k);
for ii = 1:k
    Cluster_YA{1,ii} = find(Cluster_Age{1,ii}<41);
    Cluster_MA{1,ii} = find(Cluster_Age{1,ii}>40);
end
%% Calculates Dunn Index
idx = zeros(size(Matrix,2),1);
for ii = 1:size(Cluster_Num_Cell,2)
    for a = 1:length(Cluster_Num_Cell{1,ii}{1,1})
        idx(Cluster_Num_Cell{1,ii}{1,1}(1,a),1) = ii;
    end
end

distM = zeros(size(Matrix,2),size(Matrix,2));%preallocate variable
m = 1;
while m < size(Matrix,2)+1
    for a = 1:size(Matrix,2)
        distM(m,a) = dtw(Matrix(:,m),Matrix(:,a)); 
    end
    m = m+1;
end

Dunn = dunns(10,distM,idx);

%% Total Within Sum of Squares
Within_Cluster = struct('Distance',cell(1,k));%preallocate variable
for ii = 1:k
    for a = 1:size(Cluster_Vectors_Cell{1,ii}{1,1},2)
        Within_Cluster(ii).Distance(a,1) = dtw(Cluster_Vectors_Cell{1,ii}{1,1}(:,a), Cluster_Centroid(:,ii));
    end
end

Within_Cluster_Sum = zeros(1,k);%preallocate variable
for ii = 1:k
    Within_Cluster_Sum(ii) = sum(Within_Cluster(ii).Distance);
end

Within_Cluster_Total_Sum = sum(Within_Cluster_Sum);
%% Graph clustering results
%This will be used to display the number of participants in each cluster on
%the graph 
Cluster_Total = zeros(1,k*2);
a = 1;
for ii = 1:k
    Cluster_Total(1,a) = (ii);
    a = a+1;
    Cluster_Total(1,a) = Cluster_length(1,ii);
    a = a+1;
end

%Graph the final results
if filename == 'Lumbar_Flex_Ext_Kinematics_Data_Final.xlsx'
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
hold on
for ii = 1:k
    plot(Cluster_Centroid(:,ii),'LineWidth',3.0)%plots all cluster data
end
hold on
xlim([0 100])%set x axis limit
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times New Roman','fontsize',18,'FontWeight','bold')
xline(34,'--');
xline(54,'--');
xline(73,'--');
title('K-means Clustering For Lumbar Kinematics','FontName','Times New Roman','FontSize', 28,'FontWeight','bold')%title
set(get(gca,'title'),'Position',[50 45 1.00011])%Set title position
xlabel('Trunk Flexion and Extension Cycle (%)','FontName','Times New Roman','FontSize',22,'FontWeight','bold')%xlabel name, font size, and bolded
set(get(gca,'xlabel'),'Position',[50 -131 1.00011])%set xlabel position
ylabel('Trunk Flexion/Extension Angle (°)','FontName','Times New Roman','FontSize',22,'FontWeight','bold')%ylabel name, font size, and bolded
hold on
annotation('textbox',[.11 .075 .2 .01],'String','Neutral','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.36 .075 .2 .01],'String','Max Flexion','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.53 .075 .2 .01],'String','Neutral','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.66 .075 .2 .01],'String','Hyperextension','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.885 .075 .2 .01],'String','Neutral','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
set(gcf,'color','w');%changes graph background from grey to white
legend(sprintfc('Cluster %d (n = %d)',Cluster_Total),'Location','east', 'FontSize', 18)
legend boxoff 
end

if filename == 'Pelvis_Flex_Ext_Kinematics_Data_Final.xlsx'
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
hold on
for ii = 1:k
    plot(Cluster_Centroid(:,ii),'LineWidth',3.0)%plots all cluster data
end
hold on
xlim([0 100])%set x axis limit
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times New Roman','fontsize',18,'FontWeight','bold')
xline(34,'--');
xline(54,'--');
xline(73,'--');
title('K-means Clustering For Pelvis Kinematics','FontName','Times New Roman','FontSize', 28,'FontWeight','bold')%title
set(get(gca,'title'),'Position',[50 21 1.00011])%Set title position
xlabel('Trunk Flexion and Extension Cycle (%)','FontName','Times New Roman','FontSize',22,'FontWeight','bold')%xlabel name, font size, and bolded
set(get(gca,'xlabel'),'Position',[50 -55 1.00011])%set xlabel position
ylabel('Trunk Flexion/Extension Angle (°)','FontName','Times New Roman','FontSize',22,'FontWeight','bold')%ylabel name, font size, and bolded
hold on
annotation('textbox',[.11 .075 .2 .01],'String','Neutral','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.36 .075 .2 .01],'String','Max Flexion','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.53 .075 .2 .01],'String','Neutral','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.66 .075 .2 .01],'String','Hyperextension','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
annotation('textbox',[.885 .075 .2 .01],'String','Neutral','EdgeColor','none','FontName','Times New Roman','FontSize', 14,'FontWeight','bold')
set(gcf,'color','w');%changes graph background from grey to white
legend(sprintfc('Cluster %d (n = %d)',Cluster_Total),'Location','east', 'FontSize', 18)
legend boxoff
end

figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
clr = lines(k); 
hold on
for ii = 1:k
    plot(Cluster_Vectors_Cell{1,ii}{1,1},'Color',clr(ii,:))%plots all cluster data
end
hold on
for ii = 1:k
    plot(Cluster_Centroid(:,ii),'LineWidth',3.0)%plots all cluster data
end
hold on



plot(Cluster_Vectors_Cell{1,1}{1,1})
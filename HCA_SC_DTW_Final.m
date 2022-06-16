%% Silhouette Coefficient
% calculate a = average distance of i to the points in the same cluster
% calculate b = min(average distance of i to points in another cluster)
% Silhouette coefficient = 1 - a/b 
% if a < b, (or s = b/a - 1 if a > or equal to b, not the usual case)
% typically between 0 and 1
% the closer to 1 the better

%This calculates a variable in the silhouette equation.
%This calculates the distance of each vector to other vectors within the
%same cluster.

k = Number_of_Clusters;

for ii = 1:k
    n = 1;
    while n<size(Cluster_Vectors{1,ii},2)+1
        for a = 1:size(Cluster_Vectors{1,ii},2)
            Cluster_Vector(ii).Silhouette(a,n) = dtw(Cluster_Vectors{1,ii}(:,n), Cluster_Vectors{1,ii}(:,a));
        end
        n = n+1;
    end
end

%This sorts the values in each row from smallest to greatest number.
%This is so I can delete the column with all the zeros. 
%Zeros represent the distance between the same vector, which we do not
%need.
for ii = 1:k
    for a = 1:size(Cluster_Vector(ii).Silhouette,2)
        Cluster_Vector(ii).Silhouette(a,:) = sort(Cluster_Vector(ii).Silhouette(a,:));
    end
end

%This deletes the column with the zeros
for ii = 1:k
    if length(Cluster_Vector(ii).Silhouette(:,1)) > 1
        Cluster_Vector(ii).Silhouette(:,1) = [];
    elseif length(Cluster_Vector(ii).Silhouette(:,1)) == 1
        Cluster_Vector(ii).Silhouette(:,1) = Cluster_Vector(ii).Silhouette(:,1);
    end
end

%This calculates the average of for each row vector
for ii = 1:k
    for a = 1:length(Cluster_Vector(ii).Silhouette)
        Cluster_Vector(ii).Silhouette_Ave(a,:) = mean(Cluster_Vector(ii).Silhouette(a,:));
    end
end

%This calculates the distance between each vector within each cluster to
%the centroid of the other clusters. This will help determine the nearest
%centroid each vector is closest two (other than the centroid they are already in)
for ii = 1:k
    for b = 1:k
        for a = 1:size(Cluster_Vectors{1,ii},2)
            Cluster_Dist(ii).Between(b,a) = dtw(Cluster_Vectors{1,ii}(:,a), Cluster_Ave{1,b}(:,1));
        end
    end
end

%The above section also calculated the distance between each vector and
%their own centroid which we do not need. This section will delete those
%rows because we do not need them. 
for ii = 1:k
    Cluster_Dist(ii).Between(ii,:) = [];
end

%This calculates the minimum distance values within each vector distance
%measures.
if k > 2
    for ii = 1:k
        Cluster_Dist(ii).Min = min(Cluster_Dist(ii).Between);
    end
end

%This section subtracts the minimum values just calculated with the
%Cluster_Dist.Between variable. A zero represents a vector that is closest
%to that centroid.
if k > 2
    for ii = 1:k
        for a = 1:size(Cluster_Dist(ii).Between,1)
            Cluster_Dist(ii).Sim(a,:) = Cluster_Dist(ii).Between(a,:) - Cluster_Dist(ii).Min;
        end
    end
end

%In order to have those code calculate silhouette coefficients for any
%number of centroid k, I created the Clus variable, which creates a series
%of number based on how many centroid there are.
for ii = 1:k
    Clus(ii,1:k) = 1:k;
end

%This will eliminate the centroid number of that centroid. 
%So with in the first column which is centroid 1, the 1 will be changed to
%a zero and centroid two the two will be changed to a zero and so on.if k > 2
for ii = 1:k
    for a = 1:k
        if Clus(ii,a) == ii
            Clus(ii,a) = 0;
        end
    end
end

%I want to delete the zeros within this variable so I start by sorting the
%rows
for ii = 1:k
    Clus(ii,:) = sort(Clus(ii,:));
end


%This deletes the column with the zeros
Clus(:,1) = [];

if k > 2
    for ii = 1:k
        for a = 1:size(Cluster_Dist(ii).Sim,1) 
            for b = 1:size(Cluster_Dist(ii).Sim(a,:),2)
                if Cluster_Dist(ii).Sim(a,b) == 0
                    Cluster_Dist(ii).Sim(a,b) = Clus(ii,1);
                elseif Cluster_Dist(ii).Sim(a,b) > 0
                    Cluster_Dist(ii).Sim(a,b) = Clus(ii,2);
                end
            end
        end
    end
end

%This section calculates the distance between each individual vector in one
%cluster and to each vector in the other clusters.
for ii = 1:k
    for b = 1:k
        Silhouette(b,ii).Cluster = zeros(Cluster_length(1,ii),Cluster_length(1,b));
    end
end

for ii = 1:k
    n = 1;
    while n<size(Cluster_Vectors{1,ii},2)+1
        for b = 1:k
            for a = 1:size(Cluster_Vectors{1,b},2)
                Silhouette(b,ii).Cluster(n,a) = dtw(Cluster_Vectors{1,ii}(:,n), Cluster_Vectors{1,b}(:,a));
            end
        end
        n = n+1;
    end
end

%Since I only care about the distances in vectors to the opposite cluster,
%I do not need the calculates of the distances for the same vector. This
%section deletes those sections not needed.
for ii = 1:k
    Silhouette(ii,ii).Cluster = [];
end

%This averages the distance measure calculated for Silhouette.Cluster
%variable when cluster (k) equals 2. 
if k == 2
    Silhouette_Between = [];
    for ii = 1:k
        for a = 1:size(Cluster_Vectors{1,ii},2)
            Silhouette_Between(ii).Ave(a,1) = mean(Silhouette(Clus(ii,1),ii).Cluster(a,:));
        end
    end
end

if k > 2
for ii = 1:k
    for a = 1:k
        Silhouette_Between_Matrix(a,ii) = a;
    end
end

for ii = 1:k
    Silhouette_Between_Matrix(ii,ii) = 0;
end
    
for ii = 1:k
    for b = 1:k
        if Silhouette_Between_Matrix(b,ii) ~=0
            for a = 1:size(Cluster_Vectors{1,ii},2)
            Silhouette_Between(b,ii).Ave(a,1) = mean(Silhouette(b,ii).Cluster(a,:));
            end
        end
    end
end

for ii = 1:k
    a = 1;
    for b = 1:k
        if Silhouette_Between_Matrix(b,ii) ~=0     
            SC_Between(1,ii).Total(:,a) = Silhouette_Between(b,ii).Ave;
        a = a+1;
        end
    end
end

for ii = 1:k
    for a = 1:size(Cluster_Vectors{1,ii},2)
        SC_Between(1,ii).Ave(a,1) = mean(SC_Between(1,ii).Total(a,:));
    end
end

for ii = 1:k
    for a = 1:size(Cluster_Vectors{1,ii},2)
        SC_Between(1,ii).Min(a,1) = min(SC_Between(1,ii).Total(a,:));
    end
end
end

% This calculates silhouette coefficient for each time-series based on the
% equation listed at the top of the code.
% for ii = 1:k
%     for a = 1:size(Cluster_Vectors_Cell{1,ii}{1,:},2)
%         if Cluster_Vector(ii).Silhouette_Ave(a,1)<Silhouette_Between(ii).Ave(a,1)
%             Cluster_Vector(ii).SC(a,1) = 1 - (Cluster_Vector(ii).Silhouette_Ave(a,1)/Silhouette_Between(ii).Ave(a,1));
%         elseif Cluster_Vector(ii).Silhouette_Ave(a,1)>Silhouette_Between(ii).Ave(a,1)
%             Cluster_Vector(ii).SC(a,1) = (Silhouette_Between(ii).Ave(a,1)/Cluster_Vector(ii).Silhouette_Ave(a,1))-1;
%         end
%     end
% end

if k >2
for ii = 1:k
    for a = 1:size(Cluster_Vectors{1,ii},2)
    Cluster_Vector(ii).SC(a,1) = (SC_Between(ii).Ave(a,1) - Cluster_Vector(ii).Silhouette_Ave(a,1))/max(SC_Between(ii).Ave(a,1),Cluster_Vector(ii).Silhouette_Ave(a,1));
    end
end


for ii = 1:k
    for a = 1:size(Cluster_Vectors{1,ii},2)
    Cluster_Vector(ii).SC(a,1) = (SC_Between(ii).Min(a,1) - Cluster_Vector(ii).Silhouette_Ave(a,1))/max(SC_Between(ii).Min(a,1),Cluster_Vector(ii).Silhouette_Ave(a,1));
    end
end
end

if k == 2
    for ii = 1:k
        for a = 1:size(Cluster_Vectors{1,ii},2)
            Cluster_Vector(ii).SC(a,1) = (Silhouette_Between(ii).Ave(a,1) - Cluster_Vector(ii).Silhouette_Ave(a,1))/max(Silhouette_Between(ii).Ave(a,1),Cluster_Vector(ii).Silhouette_Ave(a,1));
        end
    end
end


%This places all silhouette coefficients in one vector to be averaged
Silhouette_Total = zeros(size(Matrix,2),1);
Silhouette_Total = Cluster_Vector(1).SC;
for ii = 1:k-1
    for a = 1:length(Cluster_Vector(ii+1).SC)
    Silhouette_Total(end+1) = Cluster_Vector(ii+1).SC(a,1);
    end
end

for ii = 1:k
    Silhouette_Coefficient_Cluster{:,ii} = Cluster_Vector(ii).SC;
end

% writefile = 'Kmean_Silhouette_Coefficients_Pelvis_8Clusters.xlsx';
% for ii = 1:k
%     xlswrite(writefile,Silhouette_Coefficient_Cluster{ii},ii,'A2')
% end

%Each time-series silhouette coefficient is averaged to create an overall
%silhouette score.
Silhouette_Score = mean(Silhouette_Total);

%Calculates silhouette score for each cluster.
for ii = 1:k
    Silhouette_Score_Cluster(1,ii) = mean(Cluster_Vector(ii).SC);
end

figure(3)
barh(Silhouette_Total)
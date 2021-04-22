function [Population,W] = updateWeight4(Population,W,Z,EP,N,T)
% Delete overcrowded subproblems and add new subproblems

%% weight addition
fmax = max(EP.objs,[],1);
fmin = min(EP.objs,[],1);
if fmax == fmin
    PopObj = (Population.objs-repmat(fmin,[length(Population),1]))./repmat(fmax,[length(Population),1]);
    EPObj  = (EP.objs-repmat(fmin,[length(EP),1]))./repmat(fmax,[length(EP),1]);
else
    PopObj = (Population.objs-repmat(fmin,[length(Population),1]))./repmat(fmax-fmin,[length(Population),1]);
    EPObj  = (EP.objs-repmat(fmin,[length(EP),1]))./repmat(fmax-fmin,[length(EP),1]);
end

% calculate the radius of niche in the archive
d  = pdist2(EPObj,EPObj);
d(logical(eye(length(d)))) = inf;
sd = sort(d,2);
r  = median(sd(:,1));

index = ones(1,length(EP));
for i = 1:length(EP)
    % calculate the Euclidean distance between the Archive and population
    Dis = pdist2(EPObj(i,:),PopObj);
    if min(Dis) <= r
        index(i) = 0;
    end
end
UndevelopedS = EP(logical(index));

% according to the undeveloped q, weight gneration
if ~isempty(UndevelopedS)
    UndevelopedW = (UndevelopedS.objs - Z)./(repmat(sum(UndevelopedS.objs-Z,2),[1,size(UndevelopedS.objs,2)]));
    UndevelopedW = max(UndevelopedW,1e-6);
end

Undeveloped_Promising_S = []; 
Undeveloped_Promising_W = [];

if ~isempty(UndevelopedS)
    index = zeros(1,length(UndevelopedS));
    % find the neighbouring weight vectors and corresponding solutions
    for i = 1:size(UndevelopedW,1)
        % find the neighbouring weight vectors and corresponding solutions
        B = pdist2(UndevelopedW(i,:),W);
        [~,B] = sort(B,2);
        B = B(1:T); % neighboring weight vectors in population
        NeighS = Population(B); % neighboring solutions in population
        g_UndevelopedS = max(repmat(abs(UndevelopedS(i).obj-Z),length(B),1)./UndevelopedW(i,:),[],2);
        g_Neigh = max(abs(Population(B).objs-repmat(Z,length(B),1))./UndevelopedW(i,:),[],2);
        % compare with each solution in neighborhood
        for j = 1:T
            if g_UndevelopedS(j) < g_Neigh(j)
                index(i) = 1;
            elseif (g_UndevelopedS(j) == g_Neigh(j)) && (sum(UndevelopedS(i).obj)<sum(NeighS(j).obj))
                index(i) = 1;
            else
                index(i) = 0; 
                break;
            end
        end
    end
    Undeveloped_Promising_S = UndevelopedS(logical(index));
    Undeveloped_Promising_W = UndevelopedW(logical(index),:);
end

%% weight deletion
if ~isempty(Undeveloped_Promising_W)
    % find the number of shared weight vectors for each solution
    newPopulation = [Population,Undeveloped_Promising_S];
    newW = [W;Undeveloped_Promising_W];
    % Normalization
    PCObj = newPopulation.objs;
    fmax  = max(PCObj,[],1);
    fmin  = min(PCObj,[],1);
    PCObj = (PCObj-repmat(fmin,size(PCObj,1),1))./repmat(fmax-fmin,size(PCObj,1),1);
    
    Distance = pdist2(PCObj,PCObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,length(newPopulation));
    while sum(Del) < length(newPopulation) - N
        Remain = find(~Del);
        Temp = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
    newPopulation(Del) = [];
    newW(Del,:) = [];
    Population = newPopulation;
    W = newW;
end

end
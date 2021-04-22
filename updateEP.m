function EP = updateEP(EP,Offsprings,nEP,Z)
% Update the external population

    EP = [EP,Offsprings];
    EP = EP(NDSort(EP.objs,1)==1);
    [~,ia,~] = unique(EP.objs,'rows');
    EP = EP(ia);
    newEP=[];

    if length(EP) > nEP
        PopObj = EP.objs - repmat(Z,length(EP),1);
        [N,M] = size(PopObj);
        %% Normalization
        % Detect the extreme points
        Extreme = zeros(1,M);
        w       = zeros(M)+1e-6+eye(M);
        for i = 1 : M
            [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
        end
        % Calculate the intercepts of the hyperplane constructed by the extreme
        % points and the axes
        Hyperplane = PopObj(Extreme,:)\ones(M,1);
        a = 1./Hyperplane;
        if any(isnan(a))
            a = max(PopObj,[],1)';
        end
        % Normalization
        PopObj = PopObj./repmat(a',N,1);
        NormalizedEP_Objs = PopObj;  % normalized EP
        % 	% Normalization
        
        %-------------------determine PF shape---------q--------------%
        Cosine   = 1 - pdist2(PopObj,ones(1,size(PopObj,2)),'cosine');
        Distance = sqrt(sum(PopObj.^2,2)).*sqrt(1-Cosine.^2);
      
        center =[];
        d=0;
        for ii = 1:size(PopObj,2)
            [~,center(ii)] = min(Distance);
            Distance(center(ii))=[];
            d = d + sqrt(sum(PopObj(center(ii),:).^2,2));
        end
        q = d/sqrt(size(PopObj,2));
       
        index = zeros(1,size(PopObj,1));
        for i = 1:size(PopObj,1)
            index(i) = any(PopObj(i,:)> 1)  || any(PopObj(i,:)< 0) ;  
        end

        NormalizedS_Add_Objs = NormalizedEP_Objs(~index,:);
        S_Add = EP(~index);
 
        NormalizedS_Minus_Objs = NormalizedEP_Objs(index == 1,:); %normalized S_Minus
        S_Minus = EP(index == 1);
        
        if length(S_Add) > nEP
            newEP = [newEP,EP(Extreme)];
            Normalized_newEP_Objs = PopObj(Extreme,:);
            S_Add(ismember(NormalizedS_Add_Objs,PopObj(Extreme,:),'rows')==1) = [];
            NormalizedS_Add_Objs(ismember(NormalizedS_Add_Objs,PopObj(Extreme,:),'rows')==1,:) = [];
       
              while length(newEP) < nEP                  
                      if q > 1.1  
                      angle = acos(1-pdist2(NormalizedS_Add_Objs,Normalized_newEP_Objs,'cosine'));
                      else
                      %---------------Convex,reference point use (1,1,1,...) -------------------%
                      angle = acos(1-pdist2(NormalizedS_Add_Objs-repmat(ones(size(NormalizedS_Add_Objs,1),size(NormalizedS_Add_Objs,2)),[1,1]),Normalized_newEP_Objs-repmat(ones(size(Normalized_newEP_Objs,1),size(Normalized_newEP_Objs,2)),[1,1]),'cosine'));
                      %---------------Convex,reference point use (1,1,1,...) -------------------%
                      end
                       % maximum vector angle first
                       [~,rho] = max(min(angle,[],2));
                       newEP = [newEP,S_Add(rho)];
                       Normalized_newEP_Objs = [Normalized_newEP_Objs;NormalizedS_Add_Objs(rho,:)];
                       S_Add(rho) = [];
                       NormalizedS_Add_Objs(rho,:) = [];
                       
              end
                EP = newEP;
        else
            EP = S_Add;
        end
     
    end

end





function AMAWV(Global)
% <algorithm> <A>

% delta --- 0.9 --- The probability of choosing parents locally

%------------------------------- Reference --------------------------------
% A Novel Archive Maintenance for Adapting Weight Vectors in Decomposition
% based Multi-objective Evolutionary Algorithm.
%------------------------------- Reference --------------------------------


%% Parameter setting
[delta,~] = Global.ParameterSet(0.9,2);

%% Generate the weight vectors
[W,Global.N] = UniformlyRandomlyPoint(Global.N,Global.M);
% Size of neighborhood
T = ceil(Global.N/10);
nr = ceil(Global.N/100);
% Size of external elite
nEP = ceil(Global.N*2);

%% Detect the neighbours of each solution
B = pdist2(W,W);
[~,B] = sort(B,2);
B = B(:,1:T);

%% Generate random population
Population = Global.Initialization();
Z          = min(Population.objs,[],1);

%% Optimization
EP = updateEP(Population,[],nEP);
adaptation_moment=round(Global.maxgen*0.05);

while Global.NotTermination(Population)
    % For each solution
    Offsprings(1:Global.N) = INDIVIDUAL();
    for i = 1 : Global.N
        % Choose the parents
        if rand < delta
            P = B(i,randperm(size(B,2)));
        else
            P = randperm(Global.N);
        end
        
        % Generate an offspring
        Offsprings(i) = GAhalf(Population(P(1:2)));  % GAhalf operator
        
        % Update the ideal point
        Z = min(Z,Offsprings(i).obj);
        
        % Update the solutions in P by Tchebycheff approach
        g_old = max(abs(Population(P).objs-repmat(Z,length(P),1))./W(P,:),[],2);
        g_new = max(repmat(abs(Offsprings(i).obj-Z),length(P),1)./W(P,:),[],2);
        Population(P(find(g_old>=g_new,nr))) = Offsprings(i);
        
    end
    
    
    EP = updateEP(EP,Offsprings,nEP,Z);   % Archive maintenance
    
    
    if mod(Global.gen,adaptation_moment)==0 && Global.gen/Global.maxgen <= 0.9 && Global.gen/Global.maxgen >= 0.1
        
        if ~isempty(EP)
            [Population,W] = updateWeight4(Population,W,Z,EP,Global.N,T);
            
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);
        else
        end
    end
end
end
function MOEADBL(Global)
% <algorithm> <H-N>
% MOEA/D-BL: A Multiobjective Evolutionary Algorithm Based on Decomposition
% kind --- 1 --- The type of aggregation function

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% PSF is the best sclarizing function for this structure

%% set default operator
operator=Global.operator;
% operator=@ADE;

nParent=2;
if isequal(operator,@DE)
    nParent=3;
end

% fully vectorised cosine distance between two vectors x and y, which is 10
% times faster than pdist2 in loops
cos_dist=@(x,y) x*y'./sqrt(sum(x.^2,2)*sum(y.^2,2)');
%% Parameter setting
[kind, nr,nmlz] = Global.ParameterSet(1,2,1);
[fun,ideal_nadir]=scalarFun(kind);

%% Generate the weight vectors
% lower-level (dense) weight vectors
W=UniformPoint(Global.N,Global.M);
Global.N=size(W,1);

% upper-level (sparse) weight vectors
S=2*ceil(sqrt(Global.N));
W_g = UniformPoint(S,Global.M);

if Global.M>2
    W_g=UniformPoint(50*Global.N,Global.M);
    [~,W_g]=kmeans(W_g,S,'Distance','cosine','MaxIter',1000);
end
S=size(W_g,1);

% adjust upper-level weight vectors
dist=pdist2(W_g,W,'cosine');
[~,Io]=min(dist,[],1);
T=max(histc(Io,unique(Io)));

if  T<5
    warning('neighboorhood size is or too small, please increase upper level size');
end

for i=1:S
    Ii=find(Io==i);
    if length(Ii)<T
        try % use the mink function introduced in matlab2017
        [~,Ii]=mink(dist(i,:),T);
        catch % use the sort method if mink is not available
           [~,Ii]=sort(dist(i,:));
           Ii=Ii(1:T);  
        end
    end
    B(i,:)=Ii;
    
    % % Variant 1: no adjustment to group centres
    W_g(i,:)=mean(W(Ii,:),1); % adjust sparse weight vectors
end


% relationship between lower-level (dense) and upper-level (sparse) weight vectors
Bg = pdist2(W_g,W_g,'cosine');
[~,Bg] = sort(Bg,2);
Bg = Bg(:,1:S);


Bw = pdist2(W,W,'cosine');
[~,B_index] = sort(Bw,2);
Bl = B_index(:,1:T);

clear Del dist Io Ii B_index Bw

%% Generate random population
Population = Global.Initialization();
Z = min(Population.objs,[],1);
Z_ = max(Population.objs,[],1);

igd_curve=[];
%% Optimization
while Global.NotTermination(Population)
%     igd_curve(end+1)=mean(min(pdist2(Global.PF,Population.objs),[],2));
    alpha=Global.M^(0.5*Global.M)*max(0,Global.M-3);%(15-8*Global.evaluated/Global.evaluation);
    % For each neighborhood
    for i = 1 : S
        B1=[];s1=0;
        for k=1:S
            B1=union(B1,B(Bg(i,k),:));
            [~,ia]=uniquetol(Population(B1).objs,1e-6,'ByRows',true,'DataScale',ones(1,Global.M));
            if length(ia)>=T
                pool=B1(ia);
                pool=pool(1:T); % select the first T parents
                break;
            end
            if k==1
                s1=length(ia);
            end
        end
        % select mating pool
        for j=1 : T
            
            % Choose the parents
            %if rand < (Global.evaluated/Global.evaluation)
            %Pl=min(s1/T, 0.9);
            %if rand<0.5
               % Pl=0.9;
            %else
               % Pl=s1/T;
            %end
            
            
            if rand<(s1-1)/T
                if rand<0.5
                    P = pool(randperm(length(pool)));
                else
                    r=B(i,randi(T));
                    B1=Bl(r,:);
                    P = B1(randperm(length(B1)));
                    
                    % Variant 4: for each subproblem
                    %P = Bl(B(i,j),randperm(T));
                end
            
            % % Variant 2: use ordinary mating selection
            %if rand<0.9
            %   P = Bl(B(i,j),randperm(T)); % similar to MOEA/D-DE
            else
                P = randperm(Global.N);
            end
            
            % Generate an offspring
            Offspring = Global.Variation(Population(P(1:nParent)),1,operator);
            
            if nmlz==1
                Z_Z=Z_-Z;
            else
                Z_Z=1;
            end
            
            n_offspring=translate_obj(Offspring.objs,Z, Z_, ideal_nadir)./(Z_Z);
            %[~,best_g]=pdist2(W_g,n_offspring,'cosine','Smallest',1);
            [~,best_g]=max(cos_dist(W_g,n_offspring));
            if rand<0.5
                %[~,best_l]=pdist2(W(B(best_g,:),:),n_offspring,'cosine','Smallest',1);
                [~,best_l]=max(cos_dist(W(B(best_g,:),:),n_offspring));
                P=Bl(B(best_g,best_l),:);
                candidates=P(randperm(length(P)));
            else
                P=B(best_g,:);
                candidates=P(randperm(length(P)));
            end
            
            % Update the ideal point
            Z = min(Z,Offspring.obj);
            
            % avoid duplicate offspring
            if sum(sum((Population(candidates).objs-Offspring.objs).^2,2)<1e-10)>2*nr
                continue;
            end
            
            g_old=fun(W(candidates,:),translate_obj(Population(candidates).objs, Z, Z_, ideal_nadir)./(Z_Z), alpha);
            g_new=fun(W(candidates,:),n_offspring,alpha);
            
            
            % Update the solution
            update=find(g_old>=g_new,nr);
            
            % % Variant 3: no repair to solutions
            if isempty(update)
                if any((Offspring.objs-Z)<1.1*Z_Z)
                    % %builtin('_mergesimpts',Population(tmp).objs,[1e-3 1e-3],'first')
                    [~,~, ic]=uniquetol(Population(P).objs,min(1e-4,1/Global.N),'ByRows',true,'DataScale',ones(1,Global.M));
                    [ib,ia]=max(histc(ic,unique(ic))); % count how many times rows occur and find the max.
                    if rand<(ib-nr)/length(P)
                        ia=P(ic==ia); % subproblems with most repeated solutions
                        g_old=fun(W(ia,:),translate_obj(Population(ia).objs, Z, Z_, ideal_nadir)./(Z_Z),alpha);
                        [~,ib]=max(g_old); % find the least suitable subproblem
                        Population(ia(ib))=Offspring;
                    end
                end
            else
                Population(candidates(update)) = Offspring;
            end
        end
    end
    
    nondom_id=NDSort(Population.objs,Population.cons,1)==1;
    Z_ = max(Population(nondom_id).objs,[],1);

%     if Global.evaluated>=Global.evaluation
%         Population=Population(nondom_id);
%     end
% if T*S+Global.evaluated>Global.evaluation
% save(['MOEADBL_',func2str(Global.problem), '_run',num2str(Global.run),'_igd.mat'],'igd_curve');
% end
end
%     save(['MOEADBL_',func2str(Global.problem), '_run',num2str(Global.run),'_igd.mat'],'igd_curve');
    %save(['Data',filesep,'MOEADBL',filesep,func2str(Global.problem), '_run',num2str(Global.run),'_igd.mat'],'igd_curve');
end

function obj=translate_obj(obj,z_ideal,z_nadir, bool)
if bool % translate obj based on z_ideal
    obj=obj-z_ideal;
else
    obj=z_nadir-obj;
end
end

function [fun, ideal_nadir]=scalarFun(kind)
ideal_nadir=1;
% Select aggregating functions
switch kind
    case 1
        % PSF
        fun=@PSF;
    case 2
        % inverted PSF
        fun=@IPSF;
        ideal_nadir=0;
    case 3
        % PBI approach
        fun=@PBI;
    case 4
        % MSF
        fun=@MSF;  
    case 5
         % Tchebycheff approach
        fun=@(W, pop, beta) max(pop./W,[],2);
    case 6
        fun=@IPBI;
        ideal_nadir=0;
end
end

% Multiplicative Tchebycheff approach
function y=MSF(W,pop, beta)
tmp=pop./W;
alpha=beta*size(W,2)*min(W,[],2);
if size(alpha)==1
    y=max(tmp,[],2).^(1+alpha)./min(tmp,[],2).^(alpha);
else
    y=bsxfun(@power, max(tmp,[],2), 1+alpha)./bsxfun(@power, min(tmp,[],2), alpha);
end
end

function  y=PBI(W, pop, beta)
% PBI approach
normW   = sqrt(sum(W.^2,2));
normP   = sqrt(sum(pop.^2,2));
CosineP = sum(pop.*W,2)./normW./normP;
y   = normP.*CosineP + beta*normP.*sqrt(1-CosineP.^2);
end

function  y=IPBI(W, pop, beta)
% PBI approach
normW   = sqrt(sum(W.^2,2));
normP   = sqrt(sum(pop.^2,2));
CosineP = sum(pop.*W,2)./normW./normP;
y   = -normP.*CosineP + 1.1*normP.*sqrt(1-CosineP.^2);
end

function y=PSF(W, pop, beta)
normW   = sqrt(sum(W.^2,2));
normP   = sqrt(sum(pop.^2,2));
CosineP = sum(pop.*W,2)./normW./normP;
d1=max(pop./W,[],2);
y   = d1 + beta*normP.*sqrt(1-CosineP.^2);
end

function y=IPSF(W, pop, beta)
d1=min((pop+1e-6)./W,[],2);
y=-d1;
end

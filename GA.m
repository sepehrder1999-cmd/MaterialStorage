clc;
clear;
close all;

%% Problem Definition

global NFE;
NFE=0;

% CostFunction=@(x) MinOne(x)
model=CreateModel();

%tabe hazine tabei az f hast ba khorooji MyCost(f,model)
% CostFunction = @(sol) LayoutCost(sol,model);  

% اندازه‌ی کروموزوم (با استفاده از یک راه‌حل تصادفی)
temp = CreateRandomSolution(model);
nVar = length(temp.pos);    % Number of decision variables

% Decision Variables Matrix Size    Matrix 1 satr va nvar sotooni
VarSize = [1 nVar];   

%% GA Parameters

MaxIt = 50;      % Max No. of Iterations
nPop = 60;        % Population Size

pc = 0.9;                     % Crossover percentage
nc = 2 * round(pc * nPop / 2);      % No. of offsprings (Parents)

pm = 0.3;         % Mutation percentage
nm = round(pm * nPop);      % No. of Mutants

mu=0.04;         % Mutation Rate


ANSWER=questdlg('Choose Selection Method:','Genetic Algorithm','Roulette Wheel','Tournament','Random','Roulette Wheel');

UseRouletteWheelSelection=strcmp(ANSWER,'Roulette Wheel');
UseTournamentSelection=strcmp(ANSWER,'Tournament');
UseRandomSelection=strcmp(ANSWER,'Random');


if UseRouletteWheelSelection
    Selection_Pressure=4;         %Selection pressure
end

if UseTournamentSelection
    TournamentSize=3;   %Tournament Size
end

% pause(0.1);


%% Initialization

empty_individual.Position = [];      % Empty structure for individuals
empty_individual.Cost = [];
empty_individual.Sol=[];

pop = repmat(empty_individual, nPop, 1);

for i = 1:nPop
    % Initialize Position
    sol = CreateRandomSolution(model);
    pop(i).Position = sol.pos;
    % Evaluation
    % [pop(i).Cost ,pop(i).Sol] = CostFunction(pop(i).Position);
    disp(['pop ' num2str(i) ' has been created']);
end

% --- مرحله 1: محاسبه مقادیر خام 5S برای جمعیت اولیه ---
raw_penalties_pop = cell(nPop, 1);
for i = 1:nPop
    raw_penalties_pop{i} = calculateRaw5SPenalties(pop(i).Position, model);
end

num_stages = numel(model.Stage);
all_penalties_matrix = zeros(nPop * num_stages, 4); % 4 نوع جریمه داریم
counter = 1;
for i = 1:nPop
    penalties_per_stage = raw_penalties_pop{i}; % این یک struct array به اندازه T است
    for t = 1:num_stages
        all_penalties_matrix(counter, 1) = penalties_per_stage(t).P_Sort;
        all_penalties_matrix(counter, 2) = penalties_per_stage(t).P_SetInOrder;
        all_penalties_matrix(counter, 3) = penalties_per_stage(t).P_Shine;
        all_penalties_matrix(counter, 4) = penalties_per_stage(t).P_Standardize;
        counter = counter + 1;
    end
end

% --- مرحله 2: پیدا کردن Min و Max برای هر جریمه ---
s5_norm_params.Sort_min = min(all_penalties_matrix(:,1));
s5_norm_params.Sort_max = max(all_penalties_matrix(:,1));
s5_norm_params.SetInOrder_min = min(all_penalties_matrix(:,2));
s5_norm_params.SetInOrder_max = max(all_penalties_matrix(:,2));
s5_norm_params.Shine_min = min(all_penalties_matrix(:,3));
s5_norm_params.Shine_max = max(all_penalties_matrix(:,3));
s5_norm_params.Standardize_min = min(all_penalties_matrix(:,4));
s5_norm_params.Standardize_max = max(all_penalties_matrix(:,4));

% --- مرحله 3: ارزیابی هزینه نهایی هر فرد با پارامترهای نرمال‌سازی ---
for i = 1:nPop
    pop(i).Sol.raw_penalties = calculateRaw5SPenalties(pop(i).Position, model);
    [pop(i).Cost, pop(i).Sol] = LayoutCost(pop(i).Position, model, s5_norm_params, pop(i).Sol);
end

% Sort Population
Costs = [pop.Cost];     
% ye matrise hazine ke har sotoonesh meghdare hazine jamiato mige

[Costs, SortOrder] = sort(Costs);   
%cost moratab shode dar Costs va tartib moratab shode dar SortOrder

pop = pop(SortOrder);           
%pop be tartib chide shode

% Store best solution
BestSol = pop(1);

% Array to Hold best Cost values
BestCost = zeros(MaxIt, 1);     
% mikhaim be ezaye har iteration meghdare cost haro zakhire konim

% Store Cost
WorstCost=pop(end).Cost;

% Array to Hold Number of Function Evaluations
nfe=zeros(MaxIt,1);


%% Main Loop

figure(1);
set(gcf, 'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
for it = 1:MaxIt

    % Calculate Selection Probabilities
    if UseRouletteWheelSelection
        P=exp(-Selection_Pressure*Costs/WorstCost);
        P=P/sum(P);
    end

    % Crossover
    %vase in 2 sotoone kardim ke tedade crossover ha dastemoon biad
    popc = repmat(empty_individual, nc/2, 2); 
    for k = 1:nc/2      %tedade farzandan nesfe valedein
        if UseRandomSelection
             % entekhabe tasadofi

            i1 = randi([1 nPop]);       %yek adade random az 1 ta nPop
            i2 = randi([1 nPop]);
        end
        if UseRouletteWheelSelection
            i1=RouletteWheelSelection(P);
            i2=RouletteWheelSelection(P);
        end

        if UseTournamentSelection
            i1=TournamentSelection(pop,TournamentSize);
            i2=TournamentSelection(pop,TournamentSize);
        end

        % Select Parents
        p1=pop(i1);
        p2=pop(i2);

        % Apply Crossover
        [popc(k,1).Position, popc(k,2).Position, ~, ~] = Crossover(p1.Position, p2.Position, model);
        %position marboot be farzande aval az crossover k va posision 2vom
        %marboot be crossover k om va voroodi mogheiate 2 ta valed hast

        % % Evaluate Offsprings
        % [popc(k,1).Cost ,popc(k,1).Sol] = CostFunction(popc(k,1).Position);
        % [popc(k,2).Cost ,popc(k,2).Sol] = CostFunction(popc(k,2).Position);
    end

    popc=popc(:); % chon matris 2 sotoone bood bayad tak sotoone konim
    % popc.Position = popc.Position.pos;

    % Mutation 
    popm = repmat(empty_individual, nm, 1);
    for k = 1:nm
        % Select Parent
        i = randi([1 nPop]);
        p = pop(i);
        % ye adad random az 1 ta npop varmidare va ozve i om jamiato bar
        % midare

        % Apply Mutation
        [popm(k).Position, ~] = Mutation(p.Position, model, mu);
        
        % Evaluate Mutant
        % [popm(k).Cost ,popm(k).Sol] = CostFunction(popm(k).Position);
    end

    % Create Merged Population
    combined = [pop; popc; popm];

    % ===================================================================
    %             بخش اصلاح شده برای محاسبه و نرمال‌سازی S5
    % ===================================================================

    % --- مرحله 1: محاسبه مقادیر خام 5S برای کل جمعیت ---
    raw_penalties_pop = cell(numel(combined), 1);
    for i = 1:numel(combined)        
        raw_penalties_pop{i} = calculateRaw5SPenalties(combined(i).Position, model);
    end

    num_stages = numel(model.Stage);
    all_penalties_matrix = zeros(numel(combined) * num_stages, 4); % 4 نوع جریمه داریم
    counter = 1;
    for i = 1:nPop
        penalties_per_stage = raw_penalties_pop{i}; % این یک struct array به اندازه T است
        for t = 1:num_stages
            all_penalties_matrix(counter, 1) = penalties_per_stage(t).P_Sort;
            all_penalties_matrix(counter, 2) = penalties_per_stage(t).P_SetInOrder;
            all_penalties_matrix(counter, 3) = penalties_per_stage(t).P_Shine;
            all_penalties_matrix(counter, 4) = penalties_per_stage(t).P_Standardize;
            counter = counter + 1;
        end
    end

    
    % --- مرحله 2: پیدا کردن Min و Max برای هر جریمه ---
    s5_norm_params.Sort_min = min(all_penalties_matrix(:,1));
    s5_norm_params.Sort_max = max(all_penalties_matrix(:,1));
    s5_norm_params.SetInOrder_min = min(all_penalties_matrix(:,2));
    s5_norm_params.SetInOrder_max = max(all_penalties_matrix(:,2));
    s5_norm_params.Shine_min = min(all_penalties_matrix(:,3));
    s5_norm_params.Shine_max = max(all_penalties_matrix(:,3));
    s5_norm_params.Standardize_min = min(all_penalties_matrix(:,4));
    s5_norm_params.Standardize_max = max(all_penalties_matrix(:,4));
    % --- مرحله 3: ارزیابی هزینه نهایی هر فرد با پاس دادن پارامترهای نرمال‌سازی ---
    for i = 1:numel(combined)
        % پاس دادن پارامترهای نرمال‌سازی و مقادیر خام محاسبه شده به تابع هزینه
        combined(i).Sol.raw_penalties = calculateRaw5SPenalties(combined(i).Position, model);
        [combined(i).Cost, combined(i).Sol] = LayoutCost(combined(i).Position, model, s5_norm_params, combined(i).Sol);
    end
    
    % ===================================================================

    % Sort
    Costs = [combined.Cost];
    [Costs, SortOrder] = sort(Costs);
    combined = combined(SortOrder);
    pop = combined;

    % Truncation
    pop = pop(1:nPop); 
    %hazf mikonim migim pop mishe ozve 1 ta nPop avval
    %pop(nPop+1:end)=[] az nPop+1 ta akhare hazf mikone bejaye balaii in
    %ham mishe nevesht

    Costs=Costs(1:nPop);

   
    % pop(1).Sol = ParseSolution(pop(1).Position,model);
    % pop(1).Sol = ConstraintRepair(pop(1).Sol,model,it);
    % pop(1).Position = pop(1).Sol.pos;
    % [pop(1).Cost, pop(1).Sol] = CostFunction(pop(1).Position);
   


    % Store Best Solution ever found
    BestSol = pop(1);

    % Store Best Cost Ever Found
    BestCost(it) = BestSol.Cost;

    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);

    BestSol.Sol2 = ParseSolution(BestSol.Position, model);


    % figure(1);          % همیشه از شکل شماره 1 استفاده کن
    % clf;
    PlotSolution(BestSol.Position,model);
    drawnow;


end

%% Results

figure;
plot(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Cost');


%
% Improved Arteficial Bee Colony Algorithm in Matlab
% Warning:
% 		Any evident error in this matlab script or program idicate fault of the
% 		code author i.e. not to the algorithm
% 		
% Inspired by: 
% 		"ABC algorithm coded using C programming language" found at http://mf.erciyes.edu.tr/abc/pub/ABC.C
% References:
% 		Chun-Feng Wang, and Yong-Hong Zhang,"An Improved Artificial Bee Colony Algorithm for Solving 
%		Optimization Problems",IAENG International Journal of Computer Science, 43:3, IJCS_43_3_09,
%		27 August 2016
% @Author: Yared Tadesse.
% 			yaredtdss@gmail.com
%			July 29,2017 
%


clc
close all
clear

population = 50;
num_source = population/2;
num_parameter = 1;								
variable_size = [1 num_parameter]; 				
max_iteration = population * num_parameter*2;	
a_max = 1;
omega_min = 0.4;
omega_max = 0.9;
c1 =0.01;
c2 =0.01;
lower_bound = -10;
upper_bound = 10;
max_trial = 0.6*num_parameter*num_source;
null_bee.source = []; 
null_bee.cost = []; 
choatic_bee=repmat(null_bee,num_source,1);
bee=repmat(null_bee,num_source,1);
trial_counter = zeros(num_source,1);
source_fitness = zeros(num_source,1);
csource_fitness = zeros(num_source,1);
cost_function = @airy;%@(x)(x.^4-4*x.^2 + 2*x + 1); %@(x) sum(x.^2);%
best_food.source = 0;
best_food.cost = inf;

% CHOATIC SYSTEM INITIALIZATION
for i=1:num_source
    ch = unifrnd(0,1,1);
    ch = sin(pi * ch);
    choatic_bee(i).source = lower_bound + ch*(upper_bound - lower_bound);
    choatic_bee(i).cost = cost_function(choatic_bee(i).source);
    if choatic_bee(i).cost >= 0
        csource_fitness(i) = 1/(choatic_bee(i).cost + 1);
    elseif choatic_bee(i).cost < 0
        csource_fitness(i) = abs(choatic_bee(i).cost)+1;
    end
end

% OPPOSITION BASED LEARNING
for i=1:num_source
    bee(i).source = lower_bound + upper_bound - choatic_bee(i).source;
    bee(i).cost = cost_function(bee(i).source);
    if bee(i).cost >= 0
        source_fitness(i) = 1/(bee(i).cost + 1);
    elseif bee(i).cost < 0
        source_fitness(i) = abs(bee(i).cost)+1;
    end
    % greedy selection b/n choatic and opposition based learning soln
    if source_fitness(i) <= csource_fitness(i)
        bee(i).source = choatic_bee(i).source;
        source_fitness(i) = csource_fitness(i);
    end % endif
end % End of opposition learning

% Remembering best solution
for i=1:num_source
    if bee(i).cost < best_food.cost
        best_food.cost = bee(i).cost;
        best_food.source = bee(i).source;
    end
end
% clearing unused variable of Choatic initialization
clear choatic_bee
clear csource_fitness
clear ch
clear null_bee

% COMPLETE IABC ALGORITHM
for iteration = 1: max_iteration
    
    %EMPLOYEE BEE PHASE
    for i=1:num_source
        % declare matrix containing all number from 1 to num_source except current
        % index i
        nominee = [1:i-1 i+1:num_source];
        k = nominee(randi([1 numel(nominee)]));
        phi = a_max* unifrnd(-1,1,1); 
        % v_{ij}=x_{ij}+\phi_{ij}*(x_{ij}-x_{kj})
        mutant.source = bee(i).source + phi*(bee(i).source - bee(k).source);
        % Filtering the solution in range of boundary
        if mutant.source < lower_bound
            mutant.source = lower_bound;
        elseif mutant.source > upper_bound
            mutant.source = upper_bound;
        end
        mutant.cost = cost_function(mutant.source);
        
        % Calculating fitness function
        if mutant.cost >= 0
            mutant_fitness = 1/(mutant.cost + 1);
        elseif mutant.cost < 0
            mutant_fitness = abs(mutant.cost)+1;
        end
        
        % Apply greedy selection based on the fitness function
        if(mutant_fitness > source_fitness(i))
            bee(i).source = mutant.source;
            bee(i).cost = mutant.cost;
            source_fitness(i) = mutant_fitness;
            %reset trial_counter for the new obtained solution
            trial_counter(i) = 0;
        else
            %increment trial_counter if solution is not improved by mutant soln
            trial_counter(i) = trial_counter(i) + 1;
        end % end of greedy election
    end % End of Employee phase
    
    %calculating probablity of improvement for its selection by onLooker bee
    source_probability = source_fitness/sum(source_fitness);
    % ONLOOKER BEE PHASE
    m=1;
    % Variable n count number of succesfully choosen food sources by
    % onlooker bee
    n=1;
    while m <= num_source

        randP = unifrnd(0,1,variable_size);
       
        % Selection by on looker bee based on probability of food source
        if(randP < source_probability(n))
            % randomy select parameter to change excluding current index of food source
            nominee = [1:n-1 n+1:num_source];
            phi = a_max* unifrnd(-1,1,1); % determine random coefficient phi
            si = a_max* unifrnd(-1,1,1);
            % Selecting the global best soution based of colony
            [xtrem_value,xtrem_index] = min([bee.cost]);
            omega = omega_min + (iteration*(omega_max - omega_min))/max_iteration;
            x_subbest = best_food.source;
            x_best = bee(xtrem_index).source;
            mutant.source = omega*x_best + c1*phi*(x_best - bee(n).source) + c2*si*(x_subbest - bee(n).source);
            
            % FILTERING THE MUTANT SOLUTION
            if mutant.source < lower_bound
                mutant.source = lower_bound;
            elseif mutant.source > upper_bound
                mutant.source = upper_bound;
            end %End of extreme point filter
            
            mutant.cost = cost_function(mutant.source);
            
            % Solving for fitness of mutant solution by onlooker bee
            if mutant.cost >= 0
                mutant_fitness = 1/(mutant.cost + 1);
            elseif mutant.cost < 0
                mutant_fitness = abs(mutant.cost)+1;
            end
            
            %Apply greedy selection based on the fitness function
            if(mutant_fitness > source_fitness(n))
                bee(n).source = mutant.source;
                bee(n).cost = mutant.cost;
                source_fitness(n) = mutant_fitness;
                trial_counter(n) = 0;
            else
                %increment trial_counter if solution is not improved by mutant soln
                trial_counter(n) = trial_counter(n) + 1;
            end % End of Greedy selection for on looker bee
            m = m+1;
        end % End of probablity based selection
        n= n + 1;
        if n  > num_source
            n = 1;
        end
    end % End of Onlooker Bee
    
    % REMEMBERING BEST SOLUTION 
    for i=1:num_source
        if bee(i).cost < best_food.cost
            best_food.cost = bee(i).cost;
            best_food.source = bee(i).source;
        end
    end
    
    % SCOUT BEE Phase
    % Choatic variable ch_0 = rand(0,1) to fulfill this 
    ch_i = unifrnd(0,1,1);
    for i=1:num_source
        %Choosing the maximum number of trial iteration
        if trial_counter(i) > max_trial
            bee(i).source  = lower_bound + unifrnd(0,1)*(upper_bound - lower_bound);
        end % End of choatic based search of scout bee
    end %End of scout bee phase
end % End of allowed iteration

s=best_food.source
c=best_food.cost

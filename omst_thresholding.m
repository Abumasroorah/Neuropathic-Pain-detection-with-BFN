function [filtered_network, J_omsts] = omst_thresholding(origina_network)
[a,b] = size(origina_network);
if (a==b)
    original_network = origina_network;
else
    original_network =[origina_network;zeros(1,size(origina_network,2))];
    original_network = original_network + original_network';
end
    % Step 1: Extract OMSTs
    num_iterations = 11;  % Choose the desired number of iterations
    omsts = zeros(size(original_network, 1), size(original_network, 2), num_iterations);
    
    net_work = 1./(original_network);
    
    for i = 1:num_iterations
        % Apply Kruskal's algorithm on the inverse of the original network
        Ms_t = full(adjacency(minspantree(graph(net_work))));
        omst = Ms_t;
        omsts(:, :, i) = omst;
        
%         Ms_t2 = Ms_t + eye(size(Ms_t, 1), size(Ms_t, 2));
        Ms_t3 = 1 - Ms_t;
        
        % Substitute N-1 edges with 'Inf' in the original network
        if i < num_iterations
            net_work = net_work .* Ms_t3;
        end
    end
    
    % Step 2: Aggregate connections over OMSTs
%     aggregated_network = squeeze(sum(omsts,3));
      
    C=full(adjacency(graph(original_network)));
    Cv=sum(C(:));
%     

% Cv=2*sum(origina_network(:));
    % Step 3: Optimize the quality formula
    J_omsts = zeros(num_iterations, 1);
    for i = 1:num_iterations
        % Calculate global efficiency
%         D = distance_wei(aggregated_network);
        aggregated_network = squeeze(sum(omsts(:,:,1:i),3));
        global_efficiency(i) = efficiency_bin(aggregated_network);
        
        % Calculate cost
        cost(i) =  sum(aggregated_network(:))/Cv ;
        
        % Calculate the objective function (Global Cost Efficiency)
        J_omsts(i) = global_efficiency(i) - cost(i);
        
        % For each adding connection, update the aggregated network
%         if i < num_iterations
%             aggregated_network = aggregated_network + omsts(:, :, i);
%         end
    end
    
    % Step 4: Find the maximum value of the quality formula
    [~, max_index] = max(J_omsts);
    
    % Filter the original network based on the selected OMST
    filtered_network = omsts(:, :, 1:max_index);
    filtered_network = sum( filtered_network,3);
    filtered_network = filtered_network.*original_network;
end

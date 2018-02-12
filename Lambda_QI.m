function LambdaQI = Lambda_QI(AFs, LDs, queries, NN, error)

% AFs: allele frequencies of the individuals population, in the HapMap
% format (allele frequencies of alleles in 5th and 7th column)
% LDs: network of SNPs where the edges are weighed with the frequencies 
% of occurring together and SNPs are the nodes
% queries: # of queries asked
% NN: size of the Beacon
% error: sequencing error
% LambdaQI: the Lambda values of the queried person (as stated in the
% paper)

    LambdaQI = [];
    alreadyAsked = zeros(1, 1500);
    
    for p = 1:20
        
        markers = markersAsked(1:queries);
        MAFs = min(AFs{ismember(AFs.markerId, markers), 5}, AFs{ismember(AFs.markerId, markers), 7});
        MAFs = [AFs(ismember(AFs.markerId, markers), 3), array2table(MAFs)];
        [MAFs, ~] = sortrows(MAFs, 2);
        markers = MAFs{:, 1};
        MAFs = table2array(MAFs(:, 2));
        MAFs(MAFs==0) = 0.000000001;
        BeaconResponse = squeeze(responseALL(1:5, 1:queries));
        innersum = 0;
        
        for i = 1:queries
            if alreadyAsked(i) == 1
                continue;
            end
            if ismember(markers{i}, LDs.Nodes.Name)
                nodes = successors(LDs, markers{i});
                nodes = nodes(ismember(nodes, markers));
                MAFinner = MAFs(ismember(markers, nodes));
                nodes{end+1,1} = markers{i};
                if ~isempty(MAFinner)
                    alreadyAsked(ismember(markers, nodes(1:end-1))) = 1;
                    
                    m = size(nodes, 1);
                    summ = zeros(m, 3);
                    for j = 1:m
                        edges = findedge(LDs, nodes{j}, nodes);
                        summ(j, 1) = sum(LDs.Edges.Weight(edges(edges~=0)))/size(edges(edges~=0), 1);
                        summ(j, 2) = size(edges(edges~=0), 1);
                        summ(j, 3) = j;
                    end
                    
                    [~, b] = sort(summ(:,2));
                    summ = summ(b,:);
                    summ = summ(1:int8(m * 0.3),:);
                    [~, idx1] = max(summ(:, 1));
                    idx = summ(idx1, 3);
                    queryM = nodes{idx};
                    nodes = successors(LDs, queryM);
                    nodes = nodes(ismember(nodes, markers));                    
                    
                    MAFinner = MAFs(ismember(markers, nodes));
                    
                    if BeaconResponse(5,ismember(markersAsked(1:queries), queryM)) == 0 
                        innersum = innersum + (sum(log(((1-MAFinner(1:size(MAFinner, 1))).^2).*(error^(-1))) + log((error./(1-MAFinner(1:size(MAFinner, 1))).^2) .* (1-(1-MAFinner(1:size(MAFinner, 1))).^(2*NN))./(1-error.*(1-MAFinner(1:size(MAFinner, 1))).^(2*NN-2))).* zeros(size(MAFinner, 1), 1))*mean(LDs.Edges.Weight(findedge(LDs, queryM, nodes))));

                    else
                        innersum = innersum + (sum(log(((1-MAFinner(1:size(MAFinner, 1))).^2).*(error^(-1))) + log((error./(1-MAFinner(1:size(MAFinner, 1))).^2) .* (1-(1-MAFinner(1:size(MAFinner, 1))).^(2*NN))./(1-error.*(1-MAFinner(1:size(MAFinner, 1))).^(2*NN-2))).* ones(size(MAFinner, 1), 1))*mean(LDs.Edges.Weight(findedge(LDs, queryM, nodes))));

                    end
                end                
            end
            LambdaQI(i) = sum(log(((1-MAFs(~(squeeze(alreadyAsked(1:i))'))).^2).*(error^(-1))) + log((error./(1-MAFs(~(squeeze(alreadyAsked(1:i))'))).^2) .* (1-(1-MAFs(~(squeeze(alreadyAsked(1:i))'))).^(2*NN))./(1-error.*(1-MAFs(~(squeeze(alreadyAsked(1:i))'))).^(2*NN-2))).* BeaconResponse(5,~(squeeze(alreadyAsked(1:i))'))') + innersum;
        end
        
    end

end
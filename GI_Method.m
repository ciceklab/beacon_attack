function [answer, markersAsked] = GI_Method(Beacon, AFs, TrainingSet, CHR, h)

% Beacon: simulated beacon with SNPs as rows and individuals as columns
% CHR: victims SNP file
% TrainingSet: set of individuals that are used to train the markov chain
% AFs: allele frequencies of the individuals population, in the HapMap
% format (allele frequencies of alleles in 5th and 7th column)
% h: threshold of SNPs hidden with a MAF < h

    markersAsked = {};
    answer = [];
    order = 4;

    CHR2 = sortrows(CHR, 1);
    markers_old = CHR2{:, 1};
    snps_old = CHR2{:, 2};
    [markers, snps] = FilterBiAllelicNewMethod(AFs, markers_old, snps_old);
    snps_old(ismember(markers_old, markers)) = snps;  

    AFs_temp = AFs(AFs.referenceAlleleFrequency < h/100 | AFs.otherAlleleFrequency < h/100,:);

    MAFs = min(AFs_temp{ismember(AFs_temp.markerId, markers), 5}, AFs_temp{ismember(AFs_temp.markerId, markers), 7});
    MAFs = [AFs_temp(ismember(AFs_temp.markerId, markers), 3), array2table(MAFs)];
    [MAFs, indexMAF] = sortrows(MAFs, 2);
    MAFs = table2array(MAFs(:, 2));

    toInfer = AFs_temp(ismember(AFs_temp.markerId, markers), :);
    toInfer = toInfer(indexMAF,:);

    indexMarkersAsked = 1;
    for i = 1:size(toInfer, 1)

        index = find(ismember(TrainingSet.Properties.VariableNames, toInfer{i,3}));
        if isempty(index)
            continue;
        end
        if index <= order
            continue;
        end

        index2 = find(ismember(markers_old, toInfer{i,3}));
        if isempty(index2)
            continue;
        end
        if ismember('NN', CHR{index-order:index, 2})
            continue
        end

        train = TrainingSet(:, index-(order):index);
        train2 = high_order_file_new(train, AFs);
        if size(train2) == 1
            continue;
        end

        markovModel = probe_calculate(order, table2array(train2));
        a = permn([0 1 2], order+1);
        indiv = cell2table(CHR{index-order:index,2}');
        indiv.Properties.VariableNames = CHR{index-order:index,1};
        if sum(ismember(CHR{index-order:index-1,1}, toInfer.markerId)) > 0
            continue
        end

        indivHiO = high_order_file_new(indiv, AFs);
        cor = find(a(:,1) == indivHiO{1,1} & a(:,2) == indivHiO{2,1} & a(:,3) == indivHiO{3,1});
        if markovModel(order+1, cor) < 1
            continue
        end

        answer(indexMarkersAsked) = CheckBeacon(toInfer{i, 3}, snps_old(index2), Beacon);
        markersAsked(indexMarkersAsked) = toInfer{i, 3};
        if indexMarkersAsked >= 50
            break;
        end
        indexMarkersAsked = indexMarkersAsked + 1;
    end    
end
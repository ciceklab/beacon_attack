function [markersAsked, responseALL, Lambda] = QI_Method(Beacon, NN, AFs, CHR, queries, error)

% Beacon: simulated beacon with SNPs as rows and individuals as columns
% CHR: victims SNP file
% AFs: allele frequencies of the individuals population, in the HapMap
% format (allele frequencies of alleles in 5th and 7th column)
% queries: # of queries asked
% NN: size of the Beacon
% error: sequencing error

    Lambda = [];
    markersAsked = {};
    responseALL = [];

    CHR = sortrows(CHR, 1);
    markers = CHR{:, 1};
    snps = CHR{:, 2};
    [markers, snps] = FilterBiAllelicNewMethod(AFs, markers, snps);

    BeaconResponse(5,1) = CheckBeacon(markers(1), snps(1), Beacon);
    indexMarkersAsked = 1;
    Lambda(indexMarkersAsked) = sum(log(((1-MAFs(1:1)).^2).*(error^(-1))) + log((error./(1-MAFs(1:1)).^2) .* (1-(1-MAFs(1:1)).^(2*NN))./(1-error.*(1-MAFs(1:1)).^(2*NN-2))).* BeaconResponse(5,1:1)');
    markersAsked(indexMarkersAsked) = markers(1);
    indexMarkersAsked = indexMarkersAsked + 1;
    for i = 2:queries
        BeaconResponse(5,i) = CheckBeacon(markers(i), snps(i), Beacon);   
        Lambda(indexMarkersAsked) = sum(log(((1-MAFs(1:i)).^2).*(error^(-1))) + log((error./(1-MAFs(1:i)).^2) .* (1-(1-MAFs(1:i)).^(2*NN))./(1-error.*(1-MAFs(1:i)).^(2*NN-2))).* BeaconResponse(5,1:i)');
        markersAsked(indexMarkersAsked) = markers(i);
        indexMarkersAsked = indexMarkersAsked + 1;
    end
        responseALL(1:5, 1:queries) = BeaconResponse(:,1:queries);

end
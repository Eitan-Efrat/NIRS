extinctionFile="ExtinctionCoefficientsData.csv";
intensityData1="FN_032_V1_Postdose1_Nback.mat";
intensityData2="FN_031_V2_Postdose2_Nback.mat";
tissueFile="DPFperTissue.txt";
relDPFfile="RelativeDPFCoefficients.csv";
dx=[1 2];
SDS=3; %in cm
tissue="adult_head"; %6.26

[dHbR1, dHbO1, fig1] = CalcNIRS1('FN_031_V2_Postdose2_Nback.mat', 3, 'adult_forearm', [1 2], ...
    'ExtinctionCoefficientsData.csv', 'DPFperTissue.txt', 'RelativeDPFCoefficients.csv');

[dHbR2, dHbO2, fig2] = CalcNIRS1('FN_031_V2_Postdose2_Nback.mat', 3, 'adult_forearm', [1 2], ...
    'ExtinctionCoefficientsData.csv', 'DPFperTissue.txt', 'RelativeDPFCoefficients.csv');
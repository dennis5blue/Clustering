powerIter1 = dlmread('../../../Simulatorsrc/runSimulation/runTest/testImageScheduling/imagePower_20cam_iter1.out');
solIter1 = dlmread('../../../Simulatorsrc/runSimulation/runTest/testImageScheduling/imageSolution_20cam_iter1.out');
powerSumIter1 = sum(sum(powerIter1.*solIter1));

powerIter10 = dlmread('../../../Simulatorsrc/runSimulation/runTest/testImageScheduling/imagePower_20cam_iter10.out');
solIter10 = dlmread('../../../Simulatorsrc/runSimulation/runTest/testImageScheduling/imageSolution_20cam_iter10.out');
powerSumIter10 = sum(sum(powerIter10.*solIter10));

powerIter100 = dlmread('../../../Simulatorsrc/runSimulation/runTest/testImageScheduling/imagePower_20cam_iter100.out');
solIter100 = dlmread('../../../Simulatorsrc/runSimulation/runTest/testImageScheduling/imageSolution_20cam_iter100.out');
powerSumIter100 = sum(sum(powerIter100.*solIter100));

yaxis = [powerSumIter1 powerSumIter10 powerSumIter100];
xaxis = [1 10 100];

plot(xaxis, yaxis, 'x-.','LineWidth',2,'Color','b','DisplayName','branch-amd-bound' ); hold on;
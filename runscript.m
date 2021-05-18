mode = 1
para1 = 5
para2 = 5
concentrationM0 = modelEuler_v3(zeros(1,96), newParameters2(3546,:), randDigraphs2000{400},mode,para1,para2);
concentrationM1 = modelEulerM1_v3(zeros(1,96), newParameters2(3546,:), randDigraphs2000{400},mode,para1,para2);
concentrationM2 = modelEulerM2_v3(zeros(1,96), newParameters2(3546,:), randDigraphs2000{400},mode,para1,para2);
plotTest(concentrationM0,concentrationM1,concentrationM2, 'g1' )
plotTest(concentrationM0,concentrationM1,concentrationM2, 'g2' )
plotTest(concentrationM0,concentrationM1,concentrationM2, 'g3' )
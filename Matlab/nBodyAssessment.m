%%%%%%%%%%%%%%%%%%%%%
% N-BODY ASSESSMENT %
%%%%%%%%%%%%%%%%%%%%%
function nBodyAssessment(maxNumPointsPerNode, maxNumLevels, particleCoordinates, particleMasses)

[potential, tree]   = nbody(particleCoordinates, particleMasses, maxNumPointsPerNode, maxNumLevels);
refPotential        = bruteForce(particleCoordinates, particleCoordinates, particleMasses);

error = 100 * norm(potential - refPotential) / norm(refPotential);
fprintf('RMS percentage error = %f\n', error);

figure(2), tree.plotTree(4); axis off; hold on;
plot(particleCoordinates(1, :), particleCoordinates(2, :), 'or', 'MarkerSize', 4);

end

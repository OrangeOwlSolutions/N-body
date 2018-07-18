%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BRUTE-FORCE POTENTIAL EVALUATION FUNCTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function potential = bruteForce(targetCoordinates, sourceCoordinates, particleMasses)
% --- Brute-force N^2 potential evaluation  
  
N = size(targetCoordinates, 2);
potential = zeros(N, 1);

% --- Loop over the target points
for k = 1 : N
    
    % --- Target-source point coordinates  
    rx = targetCoordinates(1, k) - sourceCoordinates(1, :);  
    ry = targetCoordinates(2, k) - sourceCoordinates(2, :);
    r = sqrt(rx .* rx + ry .* ry);

    % --- 2D potential Green's function
    g = -1 / (2 * pi) * log(r);

    % --- Fix singularities (r = 0 -> log = inf)
    idx = find (g == inf | g == -inf);
    g(idx)=0;

    potential(k) = sum(g .* particleMasses);
end
  
end

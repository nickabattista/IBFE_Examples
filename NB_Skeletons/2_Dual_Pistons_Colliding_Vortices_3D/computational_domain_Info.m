function computational_domain_Info()

L = 2.5;                                    % length of computational domain (m)
MAX_LEVELS = 4;                           % maximum number of levels in locally refined grid
REF_RATIO  = 4;                           % refinement ratio between levels
NCOARSE = 16;                             % actual number of grid cells on coarsest grid level
N = (REF_RATIO^(MAX_LEVELS - 1))*NCOARSE; % effective number of grid cells on finest   grid level
dx = (1.0*L)/N;                           % Cartesian mesh width (m)
ds = dx*2; 

fprintf('\nDOMAIN INFO:\n',L);
fprintf('L = %f\n',L);
fprintf('Coarest N = %d\n',NCOARSE);
fprintf('Highest N = %d\n',N);
fprintf('dx = %f\n',dx);
fprintf('ds = %f\n\n',ds);

L = 1;                                    % length of computational domain (m)
MAX_LEVELS = 3;                           % maximum number of levels in locally refined grid
REF_RATIO  = 4;                           % refinement ratio between levels
NCOARSE = 32;                             % actual number of grid cells on coarsest grid level
N = (REF_RATIO^(MAX_LEVELS - 1))*NCOARSE; % effective number of grid cells on finest   grid level
dx = (1.0*L)/N;                           % Cartesian mesh width (m)
ds = dx*2; 

fprintf('\nOR DOMAIN INFO:\n',L);
fprintf('L = %f\n',L);
fprintf('Coarest N = %d\n',NCOARSE);
fprintf('Highest N = %d\n',N);
fprintf('dx = %f\n',dx);
fprintf('ds = %f\n\n',ds);


{ 

This is the first attempt to describe a split gate system with top gate

}
TITLE 'Band Structure'     
COORDINATES cartesian1  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u              { choose your own names }
SELECT  

	ERRLIM=1e-5       

DEFINITIONS    { parameter definitions }

  rho = TABLE('density_1d.dat', x)
  nx = 400

  ly = 4000
  
  well_depth = -30

  split_V = 0

  eps_0 = 1.41859713
  eps_r = 13
  eps

! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  div(eps * grad(u)) = rho { one possibility }

! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
	eps = eps_0 * eps_r
    START(0)   { Walk the domain boundary }
	POINT VALUE(u) = split_V
    LINE TO (ly) 
	POINT VALUE(u) = 0

! TIME 0 TO 1    { if time dependent }
! MONITORS         { show progress }
PLOTS            { save result displays }
	ELEVATION(rho) FROM (0) TO (ly)
	ELEVATION(u) FROM (0) TO (ly) export(nx) format '#1' file='pot.dat'
END

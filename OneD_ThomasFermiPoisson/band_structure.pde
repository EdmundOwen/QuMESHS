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

  rho = TABLE('dens.dat', x)
  nx = 500

  lx = 5000

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
    LINE TO (lx) 
	POINT NATURAL(u) = 0

! TIME 0 TO 1    { if time dependent }
! MONITORS         { show progress }
PLOTS            { save result displays }
	ELEVATION(rho) FROM (0) TO (lx)
	ELEVATION(u) FROM (0) TO (lx) export(nx) format '#1' file='pot.dat'
END

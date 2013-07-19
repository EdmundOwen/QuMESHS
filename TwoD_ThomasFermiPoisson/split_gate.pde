{
This is the first attempt to describe a split gate system with top gate
}

TITLE 'Split Gate'
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u              { choose your own names }
SELECT

	ERRLIM=1e-6

DEFINITIONS    { parameter definitions }

  rho = 0  {table('dens.dat', x, y)}
  nx = 40
  ny = 40

  lx = 1000	 ly = 4000
  lx_1 = lx+1		ly_1 = ly+1
  pmma_depth = 200
  gate_depth = 2
  top_depth = 2
  split_width = 400

  well_depth = -300
  well_width = 30

  split_V = -1
  top_V = -2

  dopent = 0.00001

  eps_0 = 1.41859713
  eps_pmma = 2.6
  eps_r = 13
  eps

! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  div(eps * grad(u))= rho { one possibility }
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
	eps = eps_0
    START(-lx/2,-ly)   { Walk the domain boundary }
	VALUE(u) = 0
	{NATURAL(u)=0}
    LINE TO (lx/2,-ly)
	NATURAL(u) = 0
	LINE TO (lx/2,pmma_depth+top_depth)
	LINE TO (-lx/2,pmma_depth+top_depth)
	NATURAL(u) = 0
	LINE TO CLOSE

  REGION 2 {Material}
	eps = eps_0 * eps_r
	START (-lx/2,-ly)
	LINE TO (-lx/2, 0)
	TO (lx/2, 0)
	TO (lx/2, -ly) TO CLOSE

  REGION 3 {PMMA Layer}
	eps = eps_0 * eps_pmma
	START (-lx/2, 0)
	LINE TO (-lx/2, pmma_depth) TO (lx/2, pmma_depth) TO (lx/2, 0) TO CLOSE

  REGION 4 {Left split gate}
	eps = eps_0
	START(-lx/2, 0)
	VALUE(u) = split_V
	LINE TO (-lx/2, gate_depth) TO (- split_width / 2, gate_depth) TO (- split_width / 2, 0) TO CLOSE
	
  REGION 5 {Right split gate}
	eps = eps_0
	START(lx/2, 0)
	VALUE(u) = split_V
	LINE TO (lx/2, gate_depth) TO (split_width / 2, gate_depth) TO (split_width / 2, 0) TO CLOSE

  REGION 6 {Top gate}
	eps = eps_0
	START(lx/2, pmma_depth)
	VALUE(u) = top_V
	LINE TO (lx/2, pmma_depth+top_depth) TO (-lx/2, pmma_depth+top_depth) TO (-lx/2, pmma_depth) TO CLOSE

  REGION 7 {Plotting region}
	eps = eps_0 * eps_r
	MESH_SPACING = 10
	START(lx/2, well_depth)
	LINE TO (lx/2, well_depth-well_width) TO (-lx/2, well_depth-well_width) TO (-lx/2, well_depth) TO CLOSE

  {REGION 8 {Substrate}
	eps = eps_0 * eps_r
	rho = dopent
	START(-lx/2, -ly)
	LINE TO (-lx/2, -ly-2000) TO (lx/2, -ly-2000) TO (lx/2, -ly) TO CLOSE}

! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
	CONTOUR(u)
PLOTS            { save result displays }
{  CONTOUR(rho)}
  CONTOUR(u)
	{CONTOUR(u) ZOOM (-lx/2, -ly, lx, ly) EXPORT(nx, ny) FORMAT '#1' FILE = 'pot.dat'}
	CONTOUR(u) ZOOM (-lx/2, well_depth-well_width, lx, well_width) {EXPORT(dens_nx, dens_ny) {FORMAT '#x#b#y#b#1'} FILE = 'pot_small.txt'}
	VECTOR(dx(u),dy(u)) ZOOM (-lx/2, -30, lx, 60)
	ELEVATION(u) FROM (-lx/2, well_depth-well_width/2) TO (lx/2,well_depth-well_width /2)
	ELEVATION(u) FROM (-100, well_depth-well_width/2) TO (100, well_depth-well_width/2) export format '#x#b#1' file = 'data.dat'
	ELEVATION(u) FROM (0, 0) TO (0, -ly)
	TRANSFER() FILE = 'mesh.dat'
END










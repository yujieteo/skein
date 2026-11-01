// The entry of the package mesher: the triangle and tetrahedron mesher of
// spec §9.2, and its quad and hex cells (grid.goose, quad.goose, sweep.goose). A
// program imports it with `import mesher from "<path to src>/mesher";`,
// builds a shape with the functions of geometry.goose, calls mesh_2d,
// mesh_3d or mesh_cells, and reads the
// mesh of mesh.goose and its mesh_quality.
namespace mesher;
import .geometry;
import .mesh;
import .quality;
import .tri;
import .tet;
import .grid;
import .refine;
import .quad;
import .sweep;
import .improve;

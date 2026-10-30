// The entry of the package mesher: the triangle and tetrahedron mesher of
// spec §9.2. A program imports it with
// `import mesher from "<path to src>/mesher";`, builds a shape with the
// functions of geometry.goose, calls mesh_2d or mesh_3d, and reads the
// mesh of mesh.goose and its mesh_quality.
namespace mesher;
import .geometry;
import .mesh;
import .quality;
import .tri;
import .tet;

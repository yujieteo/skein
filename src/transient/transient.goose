// The entry of the package transient: implicit time stepping of
// M du/dt + K u = 0, for the `transient` procedure (spec §14.2). A program
// imports it with `import transient from "<path to src>/transient";` and
// calls the functions of api.goose. The program also declares the packages
// runtime (for runtime::cg) and solvers (for the CSR operations), as
// src/transient/test.goose does.
namespace transient;
import .api;

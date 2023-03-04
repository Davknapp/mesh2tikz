
#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h> /* A collection of exemplary cmeshes */
#include <t8_forest.h>          /* forest definition and basic interface. */
#include <forest_2_tikz/forest_to_tikz.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>


int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  int                 level = 2;
  const char         *tikz_name = "t8_tikz";
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  const double        screen_width = 30.0;
  const double        screen_height = 30.0;
  const double        cam[3] = { 0.5, 0.5, 2.0 };
  const double        focus[3] = { 0.5, 0.5, 1.0 };
  const double        up[3] = { 0.0, 1.0, 0.0 };

  const double        view_width = 2.0;
  const double        view_height = 2.0;
  const double        far = 2.0;
  const int           write_sfc = 1;

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_DEBUG);

  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;

  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TET, comm, 0, 0, 0);

  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0,
                           comm);
  t8_productionf ("write tikz to %s\n", tikz_name);

  mesh2tikz (forest, tikz_name, screen_width, screen_height, cam,
                        focus, up, view_width, view_height, far, write_sfc);
  t8_forest_write_vtk (forest, "tikz_compare");
  t8_productionf ("done writing\n");
  t8_forest_unref (&forest);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}



#include <t8.h>                 /* General t8code header, always include this. */
#include <t8_cmesh.h>           /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h> /* A collection of exemplary cmeshes */
#include <t8_forest.h>          /* forest definition and basic interface. */
#include <forest_2_tikz/forest_to_tikz.hxx>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_cmesh_vtk.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>

/*
 * This file is demonstration on how to use mesh2tikz. 
 * We construct a forest with a uniform refinement, set the camera, the view-volume
 * and the screen. We also want to write the SFC of the forest. The result is written
 * into the file "t8_tikz.tikz". If you want to compare the result you can have a
 * look at the vtu-files "tikz_compare" in Paraview. 
 */


static t8_cmesh_t
construct_cmesh(sc_MPI_Comm comm){
  t8_cmesh_t cmesh;
  double quad_vertices[12] = {0, 0, 1,
                              1, 0, 1,
                              0, 1, 1,
                              1, 1, 1};
  double  tri_one_vertices[9] = {1, 1, 1,
                                 1, 0, 1,
                                 1, 1, 0};
  double  tri_two_vertices[9] = {1, 1, 1,
                                 0, 1, 1,
                                 1, 1, 0};  
  double  tri_three_vertices[9] = {0, 1, 1,
                                 0, 1, 0,
                                 1, 1, 0}; 
  t8_geometry_c      *linear_geom = t8_geometry_linear_new (3);                         
  t8_cmesh_init(&cmesh);
  t8_cmesh_set_tree_class(cmesh,0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class(cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class(cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class(cmesh, 3, T8_ECLASS_TRIANGLE);
  t8_cmesh_register_geometry(cmesh, linear_geom);
  t8_cmesh_set_tree_vertices(cmesh, 0, quad_vertices, 4);
  t8_cmesh_set_tree_vertices(cmesh, 1, tri_one_vertices, 3);
  t8_cmesh_set_tree_vertices(cmesh, 2, tri_two_vertices, 3);
  t8_cmesh_set_tree_vertices(cmesh, 3, tri_three_vertices, 3);

  t8_cmesh_commit(cmesh, comm);
  return cmesh;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm = sc_MPI_COMM_WORLD;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  int                 level = 2;
  const char         *tikz_name = "t8_tikz";
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  const double        screen_width = 30.0;
  const double        screen_height = 30.0;
  const double        cam[3] = { 1.5, 1.5, 1.5 };
  const double        focus[3] = { 1.5, 0.5, 0.5 };
  const double        up[3] = { 0.0, 0.0, 1.0 };

  const double        view_width = 3.0;
  const double        view_height = 3.0;
  const double        far = 3.0;
  const int           write_sfc = 1;
  const int           color_mpi = 1;
  int                 **mpi_colors = (int **)malloc(3*sizeof(int *));
  T8_ASSERT(mpi_colors != NULL);
  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (comm, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_DEBUG);

  for(int i = 0; i < 3; i++){
    mpi_colors[i] = (int *)malloc(3*sizeof(int));
    for(int j = 0; j < 3; j++){
      mpi_colors[i][j] = (i == j) ? 255: 0;
    }
  }

  cmesh = construct_cmesh(comm);
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0,
                           comm);

  mesh2tikz (forest, tikz_name, screen_width, screen_height, cam,
                        focus, up, view_width, view_height, far, write_sfc, color_mpi, mpi_colors);
  t8_forest_write_vtk (forest, "tikz_compare");
  
  for(int i = 2; i >= 0; i--){
    free(mpi_colors[i]);
  }
  free(mpi_colors);
  
  t8_forest_unref (&forest);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}


#include"forest_to_tikz.hxx"
#include "../matrices/transformation_matrix.hxx"
#include <t8.h>
#include <t8_element_cxx.hxx>
#include <t8_vec.h>

/**
 * Enum for clipping purpose. 
 * 
 */
typedef enum inside_view
{
  outside = 0,
  inside = 1
} inside_view;

/**
 * Check if all elements are inside the cube [0,1]^3
 * 
 */
inside_view
clipping (double point[3])
{
  for (int i = 0; i < 3; i++) {
    if (point[i] < -1 || point[i] > 1) {
      return outside;
    }
  }
  return inside;
}

int
check_write_success (FILE *file, int freturn)
{
  if (freturn > BUFSIZ) {
    if (file != NULL) {
      fclose (file);
    }
    return 0;
  }
  else {
    return 1;
  }
}

const int t8_to_tikz[2][4] = {{0,1,3,2},
                              {0,1,2, -1}};

/**
 * Write a single triangle in a tikz-file.
 * 
 * \param[in] file  The file where the triangle should be written to.
 * \param[in] coords  The coordinates of each vertex of the triangle.
 * \return    0, if writing was successfull. -1 otherwise.
 */
static int
write_2D (FILE *file, double **coords, const int num_vertex, const t8_element_shape_t  shape)
{
  int freturn = fprintf (file, "\\draw ");
  if (!check_write_success (file, freturn)) {
    t8_errorf ("Error writing 2D-element\n");
    return -1;
  }
  for(int ivertex = 0; ivertex < num_vertex; ivertex++){
    const int tikz_corner = t8_to_tikz[(int) shape - T8_ECLASS_QUAD][ivertex];
    freturn = fprintf (file, "(%3.3f, %3.3f) -- ", coords[tikz_corner][0], coords[tikz_corner][1]);
    if (!check_write_success (file, freturn)) {
      t8_errorf ("Error writing 2D-element\n");
      return -1;
    }
  }
  freturn = fprintf(file, "cycle;\n");
  if (!check_write_success (file, freturn)) {
    t8_errorf ("Error writing 2D-element\n");
    return -1;
  }
  return 0;
}

/**
 * Project an element onto the output screen and write it into a tikz-file.
 * 
 * \param[in] forest      The forest of the element
 * \param[in] ts          The tree-scheme
 * \param[in] ltree_id    The local tree id. 
 * \param[in] element     The element that we want to write.
 * \param[in] file        A pointer to a file
 * \param[in] inv_cam     The inverse camera projection matrix
 * \param[in] perspective The perspective projection matrix
 * \param[in] screen      The screen projection matrix. 
 */
void
write_element (t8_forest_t forest, t8_eclass_scheme_c *ts,
               const t8_locidx_t ltree_id, const t8_element_t *element,
               FILE *file, const double inv_cam[4][4],
               const double perspective[4][4], const double screen[3][3])
{
  const int           num_vertex = ts->t8_element_num_corners (element);
  double            **vertex_coords = T8_ALLOC (double *, num_vertex);
  /* Iterate over each vertex of the element and transform it into screen-coordinates. */
  for (int i = 0; i < num_vertex; i++) {
    vertex_coords[i] = T8_ALLOC (double, 3);
    double              coords[3];
    t8_forest_element_coordinate (forest, ltree_id, element, i, coords);
    double              cam_coord[3];
    double              perspective_coord[3];
    mat4d_vec_multi (inv_cam, coords, cam_coord);
    mat4d_vec_multi (perspective, cam_coord, perspective_coord);
    t8_mat_vec (screen, perspective_coord, 1.0, vertex_coords[i]);
  }
  /* Write the transformed triangle. */
  const t8_element_shape_t  shape = ts->t8_element_shape (element);
  if (shape == T8_ECLASS_QUAD || shape == T8_ECLASS_TRIANGLE) {
    const int           write_return = write_2D (file, vertex_coords, num_vertex, shape);
    for (int i = num_vertex - 1; i >= 0; i--) {
      T8_FREE (vertex_coords[i]);
    }
    T8_FREE (vertex_coords);
    if (write_return) {
      t8_errorf ("Error writing cell\n");
      return;
    }
  }
  else {
    t8_errorf ("Tikz-support only for Triangles");
    if (file != NULL) {
      for (int i = num_vertex - 1; i >= 0; i--) {
        T8_FREE (vertex_coords[i]);
      }
      T8_FREE (vertex_coords);
      fclose (file);
    }
    return;
  }
}

int
centroid_sfc (FILE *file, const double old_centroid[3],
              const double centroid[3], t8_locidx_t ltree_id,
              t8_locidx_t elem_id, const double screen[3][3])
{
  /* Very first element */
  int                 freturn;
  double              screen_old[3] = { 0.0 };
  double              screen_new[3] = { 0.0 };
  t8_mat_vec (screen, old_centroid, 1.0, screen_old);
  t8_mat_vec (screen, centroid, 1.0, screen_new);

  if (ltree_id == 0 && elem_id == 0) {
    return 1;
  }
  /* First element in a new tree */
  else if (elem_id == 0) {
    freturn = fprintf (file,
                       "\\draw[red, dashed] (%3.3f, %3.3f) -- (%3.3f, %3.3f);\n",
                       screen_old[0], screen_old[1],
                       screen_new[0], screen_new[1]);
  }
  else {
    freturn = fprintf (file,
                       "\\draw[red] (%3.3f, %3.3f) -- (%3.3f, %3.3f);\n",
                       screen_old[0], screen_old[1],
                       screen_new[0], screen_new[1]);
  }
  if (!check_write_success (file, freturn)) {
    t8_errorf ("Error writing triangle\n");
    return 0;
  }
  else {
    return -1;
  }
}

void
mesh2tikz (t8_forest_t forest, const char *fileprefix,
                      const double screen_width, const double screen_height,
                      const double cam[3], const double focus[3],
                      const double up[3], const double view_width,
                      const double view_height, const double far,
                      const int write_sfc)
{
  T8_ASSERT (forest != NULL);
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (fileprefix != NULL);

  const double        near = t8_vec_dist (cam, focus);

  T8_ASSERT (far > 0);
  T8_ASSERT (near > 0);
  T8_ASSERT (view_width > 0);
  T8_ASSERT (view_height > 0);
  T8_ASSERT (screen_width > 0);
  T8_ASSERT (screen_height > 0);
  T8_ASSERT (far - near > 0);
  T8_ASSERT (t8_vec_dist (cam, focus) > 0);

  FILE               *tikzfile = NULL;
  char                tikzname[BUFSIZ];

  sc_MPI_Comm comm = t8_forest_get_mpicomm (forest);
  int mpisize;
  int mpiret  = sc_MPI_Comm_size(comm, &mpisize);
  SC_CHECK_MPI(mpiret);


  /* This feature is not supposed to be executed on more than one process. */
  if (mpisize > 1) {
    t8_errorf ("Writing tikz-files is only support for serial usage.\n");
    return;
  }
  /* zero-initialize the transformation matrices */
  double              inv_cam[4][4] = { 0.0 };
  double              perspective[4][4] = { 0.0 };
  double              screen[3][3] = { 0.0 };
  t8_gloidx_t         written_cells = 0;        /* Counter for cells in the view-volume */

  /* Compute the transformation matrices */
  inverse_camera_transformation (cam, focus, up, inv_cam);
  perspective_projection (view_width, view_height, near, far, perspective);
  screen_projection (screen_width, screen_height, screen);

  /* Construct the filename. */
  int                 freturn =
    snprintf (tikzname, BUFSIZ, "%s.tikz", fileprefix);
  if (!check_write_success (tikzfile, freturn)) {
    t8_errorf ("Error constructing filename. \n");
    return;
  }
  /* Open the file. */
  tikzfile = fopen (tikzname, "w");
  if (tikzfile == NULL) {
    t8_errorf ("Unable to open file %s\n", tikzname);
    return;
  }
  else {
    t8_debugf ("Opened file %s\n", tikzname);
  }
  /* Write tikz-file header. */
  freturn = fprintf (tikzfile, "\\begin{tikzpicture}\n");
  if (!check_write_success (tikzfile, freturn)) {
    return;
  }

  double              old_centroid[3] = { 0.0 };

  /* Iterate over all trees and all elements. Check for each element if it is
   * inside the view-volume (inside the cube [-1,1]^2 after inverse camera and
   * perspective transformation.) and write the cell if true. */
  const t8_locidx_t   num_trees = t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_trees; itree++) {
    t8_eclass_scheme_c *scheme =
      t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest,
                                                                     itree));
    const t8_locidx_t   elems_in_tree =
      t8_forest_get_tree_num_elements (forest, itree);
    /* Iterate over all elements in the tree. */
    for (t8_locidx_t ielem = 0; ielem < elems_in_tree; ielem++) {
      const t8_element_t *element =
        t8_forest_get_element_in_tree (forest, itree, ielem);
      double              centroid[3];
      /* Clipping is done with the center of an element. */
      t8_forest_element_centroid (forest, itree, element, centroid);
      double              cam_centroid[3];
      double              perspective_centroid[3];
      mat4d_vec_multi (inv_cam, centroid, cam_centroid);
      mat4d_vec_multi (perspective, cam_centroid, perspective_centroid);
      if (clipping (perspective_centroid) == inside) {
        /* Actual cell writing. */
        write_element (forest, scheme, itree, element, tikzfile, inv_cam,
                       perspective, screen);
        if (write_sfc) {
          centroid_sfc (tikzfile, old_centroid, perspective_centroid, itree,
                        ielem, screen);
          old_centroid[0] = perspective_centroid[0];
          old_centroid[1] = perspective_centroid[1];
          old_centroid[2] = perspective_centroid[2];
        }
        written_cells++;
      }
    }
  }
  /* End the tikzpicture and close the file. */
  freturn = fprintf (tikzfile, "\\end{tikzpicture}\n");
  t8_productionf ("Wrote %li cells into a tikz-file\n", written_cells);
  if (!check_write_success (tikzfile, freturn)) {
    return;
  }
  if (tikzfile != NULL) {
    fclose (tikzfile);
  }
}
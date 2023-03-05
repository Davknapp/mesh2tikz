#include"forest_to_tikz.hxx"
#include "../matrices/transformation_matrix.hxx"
#include <t8.h>
#include <t8_element_cxx.hxx>
#include <t8_vec.h>
#include <t8_eclass.h>

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
 * Check if all elements are inside the cube [-1,1]^3
 * \param[in] point a 3D vector
 * \return outside if the point is outside of the cube, inside otherwise. 
 */
static inside_view
clipping (double point[3])
{
  for (int i = 0; i < 3; i++) {
    if (point[i] < -1 || point[i] > 1) {
      return outside;
    }
  }
  return inside;
}

/**
 * Check if a write operation to a char* was successfull. 
 * 
 * \param[in] file A file pointer
 * \param[in] freturn The return-value of an fprintf
 * \return int 1 if the write was successfull, 0 otherwise. 
 */
static int
check_write_success (MPI_File file, int freturn)
{
  if (freturn > BUFSIZ) {
    t8_errorf("Error filling char\n");
    if (file != NULL) {
      MPI_File_close(&file);
    }
    return 0;
  }
  else {
    return 1;
  }
}

/**
 * Check if a write operation to an MPIFILE was successfull. 
 * 
 * \param[in] file A file pointer
 * \param[in] mpiret The return-value of an mpi-operation
 * \return int 1 if the operation was successfull, 0 otherwise. 
 */
static int
check_mpi_write_success (MPI_File file, int mpiret)
{
  if (mpiret != MPI_SUCCESS) {
    t8_errorf("Error writing mpi-file\n");
    if (file != NULL) {
      MPI_File_close(&file);
    }
    return 0;
  }
  else {
    return 1;
  }
}

/**
 * A helper look-up table to translate between the t8code-vertex order and
 * the tikz-vertex order. 
 * 
 */
const int t8_to_tikz[2][4] = {{0,1,3,2},
                              {0,1,2, -1}};

/**
 * Write a single 2 dimensional element in a tikz-file.
 * 
 * \param[in] file  The file where the element should be written to.
 * \param[in] coords  The coordinates of each vertex of the 2D-element.
 * \param[in] num_vertex  The number of vertices of the element
 * \param[in] shape  The shape of the element
 * \return    0, if writing was successfull. -1 otherwise.
 */
static int
write_2D (MPI_File file, double **coords, const int num_vertex, const t8_element_shape_t  shape, const int color_mpi,
          const int mpirank)
{
  /* We start be printing \\draw and then add points until we have a full circle */
  char draw_buffer[30];
  char buffer[BUFSIZ];
  int freturn;
  if(color_mpi){
    freturn = snprintf(draw_buffer, BUFSIZ, "\\draw[fill=mpi_color_%04d!40]", mpirank);
  }
  else{
    freturn = snprintf(draw_buffer, BUFSIZ, "\\draw ", mpirank);
  }
  /* Iterate over all vertices and print the current coordinates. */
  if(num_vertex == 3){
    freturn = snprintf (buffer, BUFSIZ, "%s (%03.3f, %03.3f) -- (%03.3f, %03.3f) -- (%03.3f, %03.3f) -- cycle;\n", draw_buffer,
              coords[0][0], coords[0][1],
              coords[1][0], coords[1][1],
              coords[2][0], coords[2][1]
              ); 
  }
  else if(num_vertex == 4){
    freturn = snprintf (buffer, BUFSIZ, "%s (%03.3f, %03.3f) -- %03.3f, %03.3f) -- (%03.3f, %03.3f) -- (%03.3f, %03.3f) -- cycle;\n", draw_buffer,
              coords[0][0], coords[0][1],
              coords[1][0], coords[1][1],
              coords[2][0], coords[2][1],
              coords[3][0], coords[3][1]
              ); 
  }
  else{
    MPI_File_close(&file);
    t8_errorf("Only quadrilateral and triangular surfaces supported\n");
    return -1;
  }
  if (!check_write_success (file, freturn)) {
      t8_errorf ("Error writing 2D-element\n");
      return -1;
    }
  else{
    int mpiret = MPI_File_write_shared(file, &buffer, strlen(buffer),MPI_CHAR, MPI_STATUS_IGNORE);
    check_mpi_write_success(file, mpiret);
    return 0;
  }
}

/**
 * Write a single 2 dimensional element in a tikz-file.
 * 
 * \param[in] file  The file where the element should be written to.
 * \param[in] coords  The coordinates of each vertex of the 2D-element.
 * \param[in] num_faces  The number of vertices of the element
 * \param[in] shape  The shape of the element
 * \return    0, if writing was successfull. -1 otherwise.
 */
static int
write_3D (MPI_File file, double **coords, t8_eclass_scheme_c *ts,
          const t8_element_t *element, const int color_mpi,
          const int mpirank)
{
  const int num_faces = ts->t8_element_max_num_faces(element);
  /* Currently there are only quad and tri-faces in t8code. */
  double **face_coords = T8_ALLOC(double *, 4);
  for(int i = 0; i < 4; i++){
    face_coords[i] = T8_ALLOC_ZERO(double, 3);
  }
  /* Iterate over all faces of an element and write them. */
  for(int iface = 0; iface < num_faces; iface++){
    const t8_element_shape_t face_shape = ts->t8_element_face_shape(element, iface);
    const int num_vertices = t8_eclass_num_vertices[face_shape];
    for(int ivertex = 0; ivertex < num_vertices; ivertex++){
      /* Get the coordinates of the vertices of the current face. */
      const int elem_to_face_vertex = ts->t8_element_get_face_corner(element, iface, ivertex);
      face_coords[ivertex][0] = coords[elem_to_face_vertex][0];
      face_coords[ivertex][1] = coords[elem_to_face_vertex][1];
      face_coords[ivertex][2] = coords[elem_to_face_vertex][2];
    }
    /* Actual writing of the face. */
    int write_return = write_2D(file, face_coords, num_vertices, face_shape, color_mpi, mpirank);
    if(write_return == -1){
      for(int i = 3; i >= 0; i--){
        T8_FREE(face_coords[i]);
      }
      T8_FREE(face_coords);
      return -1;
    }
  }
  for(int i = 3; i >= 0; i--){
    T8_FREE(face_coords[i]);
  }
  T8_FREE(face_coords);
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
               MPI_File file, const double inv_cam[4][4],
               const double perspective[4][4], const double screen[3][3], 
               const int color_mpi, const int mpirank)
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
  int write_return;
  if (t8_eclass_to_dimension[shape] == 2) {
    write_return = write_2D (file, vertex_coords, num_vertex, shape, color_mpi, mpirank);
  }
  else{
    write_return = write_3D(file, vertex_coords, ts, element, color_mpi, mpirank);
  }
  for (int i = num_vertex - 1; i >= 0; i--) {
    T8_FREE (vertex_coords[i]);
  }
  T8_FREE (vertex_coords);
  if (write_return) {
    t8_errorf ("Error writing cell\n");
    return;
  }
}

/**
 * Write the SFC of the trees in the forest.
 * 
 * \param[in] file    A Filepointer to an open file
 * \param[in] old_centroid The centroid of the previous element. Arbitrary for the first element
 * \param[in] centroid The centroid of the element with id \a elem_id
 * \param[in] ltree_id The local tree-id of the element with id \a elem_id
 * \param[in] elem_id The local element id of an element in \a ltree_id
 * \param[in] screen The screen_projection matrix. 
 * \return int 0 if writing was successfull, 1, if the very first element was put in, -1 ow. 
 */
static int
centroid_sfc (MPI_File file, const double old_centroid[3],
              const double centroid[3], t8_locidx_t ltree_id,
              t8_locidx_t elem_id, const double screen[3][3])
{
  int                 freturn;
  char                buffer[BUFSIZ];
  double              screen_old[3] = { 0.0 };
  double              screen_new[3] = { 0.0 };
  /* Project the centroids onto the the screen. */
  t8_mat_vec (screen, old_centroid, 1.0, screen_old);
  t8_mat_vec (screen, centroid, 1.0, screen_new);
  if (ltree_id == 0 && elem_id == 0) {
    /* Very first element */
    return 1;
  }
  else if (elem_id == 0) {
    /* First element in a new tree. We use a dashed line to visualize a jump
     * between trees. */
    freturn = snprintf(buffer, BUFSIZ, "\\draw[red, dashed] (%3.3f, %3.3f) -- (%3.3f, %3.3f);\n",
                       screen_old[0], screen_old[1],
                       screen_new[0], screen_new[1]);
  }
  else {
    /* Draw a red line between two elements. */
    freturn = snprintf(buffer, BUFSIZ, "\\draw[red] (%3.3f, %3.3f) -- (%3.3f, %3.3f);\n",
                       screen_old[0], screen_old[1],
                       screen_new[0], screen_new[1]);
    
  }
  if (!check_write_success (file, freturn)) {
    t8_errorf ("Error writing SFC\n");
    return 0;
  }
  else {
    t8_debugf("[D] write: %i chars", strlen(buffer));
    int mpiret = MPI_File_write_shared(file, &buffer, strlen(buffer),MPI_CHAR, MPI_STATUS_IGNORE);
    check_mpi_write_success(file, mpiret);
    return -1;
  }
}

static int
construct_header(MPI_File file, const int color_mpi, int **rgb_colors, const int mpisize)
{
  char buffer[BUFSIZ]; 
  int char_return;
  int mpiret;
  if(color_mpi){
    T8_ASSERT(rgb_colors != NULL);
    for(int i = 0; i < mpisize; i++){
      char_return = snprintf(buffer, BUFSIZ, "\\definecolor{mpi_color_%04d}{RGB}{%i, %i, %i}\n",
                      i, rgb_colors[i][0], rgb_colors[i][1], rgb_colors[i][2]);
      if(!check_write_success(file, char_return)){
        t8_errorf("Error constructing header\n");
        return -1;
      }
      mpiret = MPI_File_write_shared(file, &buffer, strlen(buffer),MPI_CHAR, MPI_STATUS_IGNORE);
      if(!check_mpi_write_success(file, mpiret)){
        t8_errorf("Error constructing header\n");
      return -1;
      }
    }
  }
  char_return = snprintf(buffer, BUFSIZ, "\\pgfdeclarelayer{background}\n\\pgfdeclarelayer{foreground}\n\\pgfsetlayers{background, main, foreground}\n");
  if(!check_write_success(file, char_return)){
    t8_errorf("Error constructing header\n");
    return -1;
  }
  mpiret = MPI_File_write_shared(file, &buffer, strlen(buffer),MPI_CHAR, MPI_STATUS_IGNORE);
  if(!check_mpi_write_success(file, mpiret)){
    t8_errorf("Error constructing header\n");
  return -1;
  }

  char_return = snprintf(buffer, BUFSIZ, "\\begin{tikzpicture}\n");
  if(!check_write_success(file, char_return)){
    t8_errorf("Error constructing header\n");
    return -1;
  }
  mpiret = MPI_File_write_shared(file, &buffer, strlen(buffer),MPI_CHAR, MPI_STATUS_IGNORE);
  if(!check_mpi_write_success(file, mpiret)){
    t8_errorf("Error constructing header\n");
    return -1;
  } 
  return 0;
}

static void
layer_block(MPI_File file, const int mpirank, sc_MPI_Comm comm, char layer[BUFSIZ], int begin_end){
  MPI_Barrier(comm);
  if(mpirank == 0){
    char buffer[BUFSIZ];
    int char_return;
    if(begin_end == 0){
      char_return = snprintf(buffer, BUFSIZ, "\\begin{pgfonlayer}{%s}\n", layer);
    }
    else{ 
      char_return = snprintf(buffer, BUFSIZ, "\\end{pgfonlayer}\n");
    }
    if(!check_write_success(file, char_return)){
      t8_errorf("Error constructing header\n");
      return;
    }
    int mpiret = MPI_File_write_shared(file, &buffer, strlen(buffer),MPI_CHAR, MPI_STATUS_IGNORE);
    if(!check_mpi_write_success(file, mpiret)){
      t8_errorf("Error constructing header\n");
      return;
    } 
  }
  MPI_Barrier(comm);
}

static void
foreground_layer(t8_forest_t forest, MPI_File file,   double inv_cam[4][4],  double perspective[4][4],
                   double screen[3][3], const int write_sfc, const int color_mpi, const int mpirank, sc_MPI_Comm comm)
{
  /* Iterate over all trees and all elements. Check for each element if it is
   * inside the view-volume (inside the cube [-1,1]^2 after inverse camera and
   * perspective transformation.) and write the cell if true. */
  double              old_centroid[3] = { 0.0 };
  char foreground[11] = "foreground";
  layer_block(file, mpirank, comm, foreground, 0);
  
  
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
        if (write_sfc) {
          centroid_sfc (file, old_centroid, perspective_centroid, itree,
                        ielem, screen);
          old_centroid[0] = perspective_centroid[0];
          old_centroid[1] = perspective_centroid[1];
          old_centroid[2] = perspective_centroid[2];
        }
      }
    }
  }
  layer_block(file, mpirank, comm, foreground, 1);
}

static int
background_layer(t8_forest_t forest, MPI_File file,  double inv_cam[4][4], double perspective[4][4],
                  double screen[3][3], const int color_mpi, const int mpirank, sc_MPI_Comm comm)
{
  /* Iterate over all trees and all elements. Check for each element if it is
   * inside the view-volume (inside the cube [-1,1]^2 after inverse camera and
   * perspective transformation.) and write the cell if true. */
  char background[11] = "background";
  layer_block(file, mpirank, comm, background, 0);
  const t8_locidx_t   num_trees = t8_forest_get_num_local_trees (forest);
  t8_locidx_t written_cells = 0;
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
        write_element (forest, scheme, itree, element, file, inv_cam,
                       perspective, screen, color_mpi, mpirank);
        written_cells++;
      }
    }
  }
  layer_block(file, mpirank, comm, background, 1);
  return written_cells;
}

void
mesh2tikz (t8_forest_t forest, const char *fileprefix,
                      const double screen_width, const double screen_height,
                      const double cam[3], const double focus[3],
                      const double up[3], const double view_width,
                      const double view_height, const double far,
                      const int write_sfc, const int color_mpi, int **rgb_colors)
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

  MPI_File            tikzfile;
  MPI_Status          status;
  char                tikzname[BUFSIZ];

  sc_MPI_Comm comm = t8_forest_get_mpicomm (forest);
  int mpisize;
  int mpiret  = sc_MPI_Comm_size(comm, &mpisize);
  SC_CHECK_MPI(mpiret);

  int mpirank;
  mpiret = sc_MPI_Comm_rank(comm, &mpirank);
 
  /* zero-initialize the transformation matrices */
  double              inv_cam[4][4] = { 0.0 };
  double              perspective[4][4] = { 0.0 };
  double              screen[3][3] = { 0.0 };

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
  mpiret = MPI_File_open(comm, tikzname, MPI_MODE_WRONLY, MPI_INFO_NULL, &tikzfile);
  check_mpi_write_success(tikzfile, mpiret);


  /* Write tikz-file header. */
  if(mpirank == 0){
    const int header_return = construct_header(tikzfile, color_mpi, rgb_colors, mpisize);
  }
  /* Ensure that the file-header is written before any other process writes.*/
  
  /*Each layer handles its barries, therefore we don't need to orchestrate the process here. */
  t8_locidx_t written_cells = background_layer(forest, tikzfile, inv_cam, perspective,screen, color_mpi, mpirank, comm);
  
  foreground_layer(forest, tikzfile, inv_cam, perspective, screen, write_sfc, color_mpi, mpirank, comm);

  
  /*Ensure that the end of the file is written after all processes have finished writing. */

  MPI_Barrier(comm);
  if(mpirank == 0){
    char close_buffer[BUFSIZ] = "\\end{tikzpicture}\n";
    mpiret = MPI_File_write_shared(tikzfile, &close_buffer, strlen(close_buffer),MPI_CHAR, MPI_STATUS_IGNORE);
    check_mpi_write_success(tikzfile, mpiret);
  }
  /*Ensure that the end of the file is written before the file is closed. */
  MPI_Barrier(comm);
  t8_productionf ("Wrote %li cells into a tikz-file\n", written_cells);
  /* Clost the file*/
  if (tikzfile != NULL) {
    MPI_File_close (&tikzfile);
  }
  
}
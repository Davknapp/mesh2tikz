#ifndef FOREST_TO_TIKZ
#define FOREST_TO_TIKZ

#include <t8_forest.h>

 /**
  * Defining a view-volume, a camera-position and an output-sreen we project 
  * the elements in the view-volume onto the screen using a pin-hole camera-modell and
  * write the 2D into a tikz-file that can be used in latex-documents.
  * 
  * This feature is not meant to run on a highly parallel computation, it is meant to
  * visualize concepts to explain details about our code. Therefore it can not run on 
  * more than a single process. Even though it is possible, we do not recomment to use
  * it for very large forests either.
  * 
  * If you want to visualize your massive parralel runs on billions of elements, we 
  * recommend to use our vtk-output, see \a t8_forest_vtk.h for more information.
  * 
  * \param[in] forest           A commited forest.
  * \param[in] fileprefix       A filename, where the tikz-file should be created
  * \param[in] screen_width     The width of the output-screen.
  * \param[in] screen_height    The height of the output-screen.
  * \param[in] cam              The position of the camera.
  * \param[in] focus            The point the camera is facing.
  * \param[in] up               The up-direction of the camera. 
  * \param[in] view_width       The width of the frontal view-plane orthogonal to the view-direction, defining the view-volume.
  * \param[in] view_height      The height of the frontal view-plane orthogonal to the view-direction, defining the view-volume.
  * \param[in] far              The distance between the camera and the back of the view-volume. 
  */

void                mesh2tikz (t8_forest_t forest,
                                const char *fileprefix,
                                const double screen_width,
                                const double screen_height,
                                const double cam[3],
                                const double focus[3],
                                const double up[3],
                                const double view_width,
                                const double view_height,
                                const double far,
                                const int write_sfc);
#endif /* FOREST_TO_TIKZ */

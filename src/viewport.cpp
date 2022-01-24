#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan; 

  Matrix3x3 translationMat1 = Matrix3x3::identity();
  translationMat1(0, 2) = -this->centerX;
  translationMat1(1, 2) = -this->centerY;
  Matrix3x3 scaleMat = Matrix3x3::identity();
  scaleMat(0, 0) = 0.5 / this->vspan;
  scaleMat(1, 1) = 0.5 / this->vspan;
  Matrix3x3 translationMat2 = Matrix3x3::identity();
  translationMat2(0, 2) = 0.5;
  translationMat2(1, 2) = 0.5;
  svg_2_norm = translationMat2 * scaleMat * translationMat1;
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462

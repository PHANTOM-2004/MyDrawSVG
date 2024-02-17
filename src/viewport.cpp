#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox(float centerX, float centerY, float vspan) {

  // Task 5 (part 2):
  // Set svg coordinate to normalized device coordinate transformation. Your
  // input arguments are defined as normalized SVG canvas coordinates.
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan;
  //这里重载了逗号，通过逗号获得矩阵的元素
  Matrix3x3 M = Matrix3x3::identity();
  M(0, 0) = 0.5 / vspan;
  M(1, 1) = 0.5 / vspan;
  //M(2, 2) = 1;
  M(0, 2) = 0.5 - 0.5 * centerX / vspan;
  M(1, 2) = 0.5 - 0.5 * centerY / vspan;

  set_svg_2_norm(M);
}

void ViewportImp::update_viewbox(float dx, float dy, float scale) {

  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;
  set_viewbox(centerX, centerY, vspan);
}

} // namespace CMU462

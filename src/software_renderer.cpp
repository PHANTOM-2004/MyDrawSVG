#include "software_renderer.h"

#ifdef __MYIMP__
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "triangulation.h"

using namespace std;

namespace CMU462 {

// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg(SVG &svg) {
  // set top level transformation
  transformation = svg_2_screen;
  //看到这里应该能明白，最后一层是svg_2_screen，实际上应该是canvas -> screen
  //的tranform，
  //(screen<-canvas)*(canvas<-local coordinate)*(local coordinate<-child coordinate)

  // draw all elements
  for (size_t i = 0; i < svg.elements.size(); ++i) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(0, 0));
  a.x--;
  a.y--;
  Vector2D b = transform(Vector2D(svg.width, 0));
  b.x++;
  b.y--;
  Vector2D c = transform(Vector2D(0, svg.height));
  c.x--;
  c.y++;
  Vector2D d = transform(Vector2D(svg.width, svg.height));
  d.x++;
  d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();
}

void SoftwareRendererImp::set_sample_rate(size_t sample_rate) {
  // Task 4:
  // You may want to modify this for supersample support
  this->sample_rate = sample_rate;
  this->supersample = sample_rate > 1;
  printf("set sample rate %zu\n", sample_rate);
  if (supersample) {
    supersample_target.resize(
        4 * target_w * target_h * sample_rate * sample_rate, 0);
  } else
    supersample_target.clear();
}

void SoftwareRendererImp::set_render_target(unsigned char *render_target,
                                            size_t width, size_t height) {
  // Task 4:
  // You may want to modify this for supersample support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;
  printf("set target %zu %zu\n", width, height);
  if (supersample) {
    supersample_target.resize(
        4 * target_w * target_h * sample_rate * sample_rate, 0);
  } else
    supersample_target.clear();
}

void SoftwareRendererImp::draw_element(SVGElement *element) {
  // Task 5 (part 1):
  // Modify this to implement the transformation stack

  //这里需要让tranformation矩阵恢复到draw之前的样子
  const auto old_tranformation = this->transformation;
  this->transformation =   this->transformation * element->transform ;

  switch (element->type) {
  case POINT:
    draw_point(static_cast<Point &>(*element));
    break;
  case LINE:
    draw_line(static_cast<Line &>(*element));
    break;
  case POLYLINE:
    draw_polyline(static_cast<Polyline &>(*element));
    break;
  case RECT:
    draw_rect(static_cast<Rect &>(*element));
    break;
  case POLYGON:
    draw_polygon(static_cast<Polygon &>(*element));
    break;
  case ELLIPSE:
    draw_ellipse(static_cast<Ellipse &>(*element));
    break;
  case IMAGE:
    draw_image(static_cast<Image &>(*element));
    break;
  case GROUP:
    draw_group(static_cast<Group &>(*element));
    break;
  default:
    break;
  }
  this->transformation = old_tranformation;
}

// Primitive Drawing //

void SoftwareRendererImp::draw_point(Point &point) {
  Vector2D p = transform(point.position);
  rasterize_point(p.x, p.y, point.style.fillColor);
}

void SoftwareRendererImp::draw_line(Line &line) {
  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line(p0.x, p0.y, p1.x, p1.y, line.style.strokeColor);
}

void SoftwareRendererImp::draw_polyline(Polyline &polyline) {
  Color c = polyline.style.strokeColor;

  if (c.a != 0) {
    int nPoints = polyline.points.size();
    for (int i = 0; i < nPoints - 1; i++) {
      Vector2D p0 = transform(polyline.points[(i + 0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i + 1) % nPoints]);
      rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
    }
  }
}

void SoftwareRendererImp::draw_rect(Rect &rect) {
  Color c;

  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(x, y));
  Vector2D p1 = transform(Vector2D(x + w, y));
  Vector2D p2 = transform(Vector2D(x, y + h));
  Vector2D p3 = transform(Vector2D(x + w, y + h));

  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0) {
    rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
    rasterize_triangle(p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c);
  }

  // draw outline
  c = rect.style.strokeColor;
  if (c.a != 0) {
    rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
    rasterize_line(p1.x, p1.y, p3.x, p3.y, c);
    rasterize_line(p3.x, p3.y, p2.x, p2.y, c);
    rasterize_line(p2.x, p2.y, p0.x, p0.y, c);
  }
}

void SoftwareRendererImp::draw_polygon(Polygon &polygon) {
  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if (c.a != 0) {
    // triangulate
    vector<Vector2D> triangles;
    triangulate(polygon, triangles);

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if (c.a != 0) {
    int nPoints = polygon.points.size();
    for (int i = 0; i < nPoints; i++) {
      Vector2D p0 = transform(polygon.points[(i + 0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i + 1) % nPoints]);
      rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
    }
  }
}

void SoftwareRendererImp::draw_ellipse(Ellipse &ellipse) {
  // Extra credit
}

void SoftwareRendererImp::draw_image(Image &image) {
  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image(p0.x, p0.y, p1.x, p1.y, image.tex);
}

void SoftwareRendererImp::draw_group(Group &group) {
  for (size_t i = 0; i < group.elements.size(); ++i) {
    draw_element(group.elements[i]);
  }
}

// Rasterization //
void SoftwareRenderer::fill_rgb(uint8_t *target, const Color &color) {
  target[0] = (uint8_t)(0xff * color.r);
  target[1] = (uint8_t)(0xff * color.g);
  target[2] = (uint8_t)(0xff * color.b);
  target[3] = (uint8_t)(0xff * color.a);
}

void SoftwareRenderer::blend(uint8_t *dst, const Color &c) {
  dst[0] = (uint8_t)(c.r + (1 - c.a) * c.r) * 255;
  dst[1] = (uint8_t)(c.g + (1 - c.a) * c.g) * 255;
  dst[2] = (uint8_t)(c.b + (1 - c.a) * c.b) * 255;
  dst[3] = (uint8_t)(c.a + (1 - c.a) * c.a) * 255;
}
// The input arguments in the rasterization functions
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point(float x, float y, Color color, const bool
                                          basic_render) {
  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= target_w)
    return;
  if (sy < 0 || sy >= target_h)
    return;

  // fill sample - NOT doing alpha blending!
  /*
  render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (color.r * 255);
  render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (color.g * 255);
  render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (color.b * 255);
  render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (color.a * 255);
  */
  if (!supersample) {
    const auto target = 4 * (sx + sy * target_w);
    fill_rgb(target + render_target, color);
    // blend(target + render_target, color);
    return;
  }else if(basic_render){
    //看似放到这里和三角形里面没有什么区别， 但是在函数内部实际上经过了边界检查，这也就是为什么在
    //三角形内部会出错
    sx = (int) floor(x * sample_rate);
    sy = (int) floor(y * sample_rate);
    fill_rgb(&supersample_target[4 * (sx + sy * target_w * sample_rate)], color);
    return;
  }

  sx *= sample_rate, sy *= sample_rate;
  for (int i = 0; i < sample_rate; i++) {
    for (int j = 0; j < sample_rate; j++) {
      const auto target = (4 * ((sx + i) + (sy + j) * target_w * sample_rate));
      fill_rgb(&supersample_target[target], color);
      // blend(&supersample_target[target], color);
    }
  }
}

template <typename VST>
static void goAlongline(float x0, float y0, float x1, float y1, VST &function) {
  const float slope = (y0 - y1) / (x0 - x1);
  // 这两步细节也很重要，顺序不能颠倒，并且需要注意坐标是如何选择的。
  // 原本是采用浮点表示，这里进行优化，使用整数点，以提高效率
  // 不用size_t就不要用size_t
  if (std::abs(slope) > 1) {
    std::swap(x0, y0);
    std::swap(x1, y1);
  }
  if (x0 > x1 && slope >= 0 || x0 < x1 && slope < 0) {
    std::swap(x0, x1);
    std::swap(y0, y1);
  }

  const float m = std::abs((y1 - y0) / (x1 - x0));
  auto x = (int)std::round(x0), y = (int)std::round(y0);
  float eps = 0.0f;

  if (slope > 1) {
    for (; x <= (int)x1; x++) {
      function(y, x);
      if (m + eps < 0.5)
        eps = eps + m; // 否则还是选择y
      else
        y++, eps = eps + m - 1;
    }

  } else if (slope >= 0) {
    for (; x <= (int)x1; x++) {
      function(x, y);
      if (m + eps < 0.5)
        eps = eps + m; // 否则还是选择y
      else
        y++, eps = eps + m - 1;
    }

  } else if (slope < -1) {
    for (; x >= (int)x1; x--) {
      function(y, x);
      if (m + eps < 0.5)
        eps = eps + m; // 否则还是选择y
      else
        y++, eps = eps + m - 1;
    }

  } else {
    for (; x > (int)x1; x--) {
      function(x, y);
      if (m + eps < 0.5)
        eps = eps + m; // 否则还是选择y
      else
        y++, eps = eps + m - 1;
    }
  }
}

void SoftwareRendererImp::rasterize_line(float x0, float y0, float x1, float y1,
                                         Color color) {
  // Task 2:
  // Implement line rasterization
#if 1
  const float sx0 = std::floor(x0), sy0 = std::floor(y0);
  const float sx1 = std::floor(x1), sy1 = std::floor(y1);
  if (sx0 == sx1 && sy0 == sy1) {
    return;
  }
#endif

  // Bresenham's algorithm in drawing line.
  //
  // printf("line\t");
  auto vst = [this, &color](float x, float y) { rasterize_point(x, y, color); };
  goAlongline(x0, y0, x1, y1, vst);
}

static constexpr float EPS = 1e-5f;
/*@brief 返回右侧的x坐标
 * @param leftx 左侧的y坐标
 */

static bool is_triangle(float x0, float y0, float x1, float y1, float x2,
                        float y2) {
  const auto squ = [](const float d) -> float { return d * d; };
  const float l1 = std::sqrt(squ(x0 - x1) + squ(y0 - y1));
  const float l2 = std::sqrt(squ(x0 - x2) + squ(y0 - y2));
  const float l3 = std::sqrt(squ(x2 - x1) + squ(y2 - y1));
  return l1 + l2 > l3 && l1 + l3 > l2 && l2 + l3 > l1;
}

static bool point_in_triangle(float px, float py, float x0, float y0, float x1,
                              float y1, float x2, float y2) {
  const Vector2D AB{x1 - x0, y1 - y0}, BC{x2 - x1, y2 - y1},
      CA{x0 - x2, y0 - y2};
  const Vector2D AP{px - x0, py - y0}, BP{px - x1, py - y1},
      CP{px - x2, py - y2};
  const float A{(float)cross(AB, AP)}, B{(float)cross(BC, BP)},
      C{(float)cross(CA, CP)};
  const bool sa = A >= 0, sb = B >= 0, sc = C >= 0;
  return sa == sb && sb == sc;
}

void SoftwareRendererImp::rasterize_triangle(float x0, float y0, float x1,
                                             float y1, float x2, float y2,
                                             Color color) {
  // Task 3:
  // Implement triangle rasterization
  // 边界情况，不构成三角形
  // printf("enter triangle \n");
  if (!is_triangle(x0, y0, x1, y1, x2, y2))
    return;

  //这里不能用size_t啊。。
  const auto x_min = (int)std::min({x0, x1, x2});
  const auto y_min = (int)std::min({y0, y1, y2});
  const auto x_max = (int)std::max({x0, x1, x2});
  const auto y_max = (int)std::max({y0, y1, y2});

  if (!supersample) {
    for (int x = x_min; x <= x_max; x++)
      for (int y = y_min; y <= y_max; y++) {
        if (point_in_triangle(x + 0.5f, y + 0.5f, x0, y0, x1, y1, x2, y2))
          rasterize_point(x + 0.5f, y + 0.5f, color);
      }
    // printf("simple triangle\n");
    return;
  }

  //printf("supersample triangle\n");
  const float box_w = 1.0/sample_rate;
  for (int x = x_min; x <= x_max; x++)
    for (int y = y_min; y <= y_max; y++) {

      for (int i = 0; i < sample_rate; i++)
        for (int j = 0; j < sample_rate; j++) {
          const float tx = x + (i + 0.5f) * box_w;
          const float ty = y + (j + 0.5f) * box_w;

          if (point_in_triangle(tx, ty, x0, y0, x1, y1, x2, y2))
            rasterize_point(tx, ty, color, true);
        }
    }

}

void SoftwareRendererImp::rasterize_image(float x0, float y0, float x1,
                                          float y1, Texture &tex) {
  // Task 6:
  // Implement image rasterization
}

// resolve samples to render target
void SoftwareRendererImp::resolve(void) {
  // Task 4:
  // Implement supersample
  // You may also need to modify other functions marked with "Task 4".
  if (!supersample)
    return;

  //printf("enter resolve\t render_target %p\t render_target_size %zu\n",
         //render_target, target_h * target_w);
  //printf("%zu \t %zu\n", target_w, target_h);
  // memset(render_target, 0, target_h * target_w);
  // size_t tot = 0;
  const int N = sample_rate * sample_rate;
  for (size_t x = 0; x < target_w; x++)
    for (size_t y = 0; y < target_h; y++) {
      int sum_a{}, sum_r{}, sum_g{}, sum_b{};
      for (size_t p = 0; p < sample_rate; p++) {
        for (size_t q = 0; q < sample_rate; q++) {
          // tot++;
          const size_t target_supersample =
              4 * ((x * sample_rate + p) +
                   (y * sample_rate + q) * target_w * sample_rate);
          sum_r += supersample_target[target_supersample + 0];
          sum_g += supersample_target[target_supersample + 1];
          sum_b += supersample_target[target_supersample + 2];
          sum_a += supersample_target[target_supersample + 3];
        }
      }

      sum_r /= N, sum_g /= N, sum_b /= N, sum_a /= N;
      const size_t r_target = 4 * (x + y * target_w);
      // printf("r_target : \t%zu \t%zu\n", r_target / 4, r_target);
      // if (r_target > target_h * target_w * 4) printf("error\t");
      render_target[r_target + 0] = sum_r;
      render_target[r_target + 1] = sum_g;
      render_target[r_target + 2] = sum_b;
      render_target[r_target + 3] = sum_a;
    }
  // printf("total resolve : %zu\n", tot);
}

} // namespace CMU462
#endif
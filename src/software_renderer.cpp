#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462
{

  // Implements SoftwareRenderer //

  void SoftwareRendererImp::draw_svg(SVG &svg)
  {

    if(sample_target) {
      memset(sample_target, 0, 4 * target_w * target_h * sample_rate * sample_rate);
    } else {
      memset(render_target, 0, 4 * target_w * target_h);
    }
    // set top level transformation
    transformation = svg_2_screen;
    while(!transformation_stack.empty()) {
      transformation_stack.pop();
    }
    // draw all elements
    for (size_t i = 0; i < svg.elements.size(); ++i)
    {
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

  void SoftwareRendererImp::set_sample_rate(size_t sample_rate)
  {

    // Task 4:
    // You may want to modify this for supersampling support
    if(this->sample_rate == sample_rate) {
      return;
    }
    this->sample_rate = sample_rate;
    if(sample_target) {
      delete [] sample_target;
      sample_target = nullptr;
    }
    if (sample_rate != 1) {
      sample_target = new unsigned char[4 * target_w * target_h * sample_rate * sample_rate];
    }
  }

  void SoftwareRendererImp::set_render_target(unsigned char *render_target,
                                              size_t width, size_t height)
  {

    // Task 4:
    // You may want to modify this for supersampling support
    this->render_target = render_target;
    if(sample_rate != 1 && (this->target_w != width || this->target_h != height)) {
      delete [] sample_target;
      sample_target = new unsigned char[4 * width * height * sample_rate * sample_rate];
    }
    this->target_w = width;
    this->target_h = height;

  }

  void SoftwareRendererImp::draw_element(SVGElement *element)
  {

    // Task 5 (part 1):
    // Modify this to implement the transformation stack

    transformation_stack.push(transformation);
    transformation = transformation * element->transform;
    switch (element->type)
    {
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
    transformation = transformation_stack.top();
    transformation_stack.pop();
  }

  // Primitive Drawing //

  void SoftwareRendererImp::draw_point(Point &point)
  {

    Vector2D p = transform(point.position);
    rasterize_point(p.x, p.y, point.style.fillColor);
  }

  void SoftwareRendererImp::draw_line(Line &line)
  {

    Vector2D p0 = transform(line.from);
    Vector2D p1 = transform(line.to);
    rasterize_line(p0.x, p0.y, p1.x, p1.y, line.style.strokeColor);
  }

  void SoftwareRendererImp::draw_polyline(Polyline &polyline)
  {

    Color c = polyline.style.strokeColor;

    if (c.a != 0)
    {
      int nPoints = polyline.points.size();
      for (int i = 0; i < nPoints - 1; i++)
      {
        Vector2D p0 = transform(polyline.points[(i + 0) % nPoints]);
        Vector2D p1 = transform(polyline.points[(i + 1) % nPoints]);
        rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
      }
    }
  }

  void SoftwareRendererImp::draw_rect(Rect &rect)
  {

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
    if (c.a != 0)
    {
      rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
      rasterize_triangle(p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c);
    }

    // draw outline
    c = rect.style.strokeColor;
    if (c.a != 0)
    {
      rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
      rasterize_line(p1.x, p1.y, p3.x, p3.y, c);
      rasterize_line(p3.x, p3.y, p2.x, p2.y, c);
      rasterize_line(p2.x, p2.y, p0.x, p0.y, c);
    }
  }

  void SoftwareRendererImp::draw_polygon(Polygon &polygon)
  {

    Color c;

    // draw fill
    c = polygon.style.fillColor;
    if (c.a != 0)
    {

      // triangulate
      vector<Vector2D> triangles;
      triangulate(polygon, triangles);

      // draw as triangles
      for (size_t i = 0; i < triangles.size(); i += 3)
      {
        Vector2D p0 = transform(triangles[i + 0]);
        Vector2D p1 = transform(triangles[i + 1]);
        Vector2D p2 = transform(triangles[i + 2]);
        rasterize_triangle(p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c);
      }
    }

    // draw outline
    c = polygon.style.strokeColor;
    if (c.a != 0)
    {
      int nPoints = polygon.points.size();
      for (int i = 0; i < nPoints; i++)
      {
        Vector2D p0 = transform(polygon.points[(i + 0) % nPoints]);
        Vector2D p1 = transform(polygon.points[(i + 1) % nPoints]);
        rasterize_line(p0.x, p0.y, p1.x, p1.y, c);
      }
    }
  }

  void SoftwareRendererImp::draw_ellipse(Ellipse &ellipse)
  {

    // Extra credit
  }

  void SoftwareRendererImp::draw_image(Image &image)
  {

    Vector2D p0 = transform(image.position);
    Vector2D p1 = transform(image.position + image.dimension);

    rasterize_image(p0.x, p0.y, p1.x, p1.y, image.tex);
  }

  void SoftwareRendererImp::draw_group(Group &group)
  {

    for (size_t i = 0; i < group.elements.size(); ++i)
    {
      draw_element(group.elements[i]);
    }
  }

  // Rasterization //

  // The input arguments in the rasterization functions
  // below are all defined in screen space coordinates

  void SoftwareRendererImp::rasterize_point(float x, float y, Color color)
  {

    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= target_w)
      return;
    if (sy < 0 || sy >= target_h)
      return;

    unsigned char* target;
    if(sample_rate != 1) {
      target = sample_target;
    } else {
      target = render_target;
    }
    // fill sample - NOT doing alpha blending!
    fill_pixel(target, sx, sy, color);
  }

  void SoftwareRendererImp::rasterize_line(float x0, float y0,
                                           float x1, float y1,
                                           Color color)
  {

    // Task 2:
    // Implement line rasterization
    unsigned char* target;
    if(sample_rate != 1) {
      target = sample_target;
    } else {
      target = render_target;
    }

    float t1 = -x0 / (x1 - x0);
    float t2 = (target_w - x0) / (x1 - x0);
    if (t1 > t2)
    {
      swap(t1, t2);
    }

    float t3 = -y0 / (y1 - y0);
    float t4 = (target_h - y0) / (y1 - y0);
    if (t3 > t4)
    {
      swap(t3, t4);
    }

    float tmin = fmax(t1, t3);
    tmin = fmax(0.0, tmin);
    float tmax = fmin(t2, t4);
    tmax = fmin(tmax, 1.0);
    if (tmin > tmax)
    {
      return;
    }

    float x0_new = x0 * (1.0 - tmin) + x1 * tmin;
    float y0_new = y0 * (1.0 - tmin) + y1 * tmin;
    float x1_new = x0 * (1.0 - tmax) + x1 * tmax;
    float y1_new = y0 * (1.0 - tmax) + y1 * tmax;

    int sx0 = (int)floor(x0_new);
    int sy0 = (int)floor(y0_new);
    int sx1 = (int)floor(x1_new);
    int sy1 = (int)floor(y1_new);

    int dx = sx1 - sx0;
    int dy = sy1 - sy0;

    if (abs(dx) >= abs(dy))
    {
      if (dx < 0)
      {
        swap(sx0, sx1);
        swap(sy0, sy1);
        dx = -dx;
        dy = -dy;
      }
      float esp = 0.0;
      int y = sy0;
      float m = (float)dy / dx;
      if (dy >= 0)
      {
        for (int x = sx0; x <= sx1; ++x)
        {
          fill_pixel(target, x, y, color);
          esp += m;
          if (esp > 0.5)
          {
            y++;
            esp -= 1.0;
          }
        }
      }
      else
      {
        for (int x = sx0; x <= sx1; ++x)
        {
          fill_pixel(target, x, y, color);
          esp += m;
          if (esp < -0.5)
          {
            y--;
            esp += 1.0;
          }
        }
      }
    }
    else
    {
      if (dy < 0)
      {
        swap(sx0, sx1);
        swap(sy0, sy1);
        dx = -dx;
        dy = -dy;
      }
      float esp = 0.0;
      int x = sx0;
      float m = (float)dx / dy;
      if (dx >= 0)
      {
        for (int y = sy0; y <= sy1; ++y)
        {
          fill_pixel(target, x, y, color);
          esp += m;
          if (esp > 0.5)
          {
            x++;
            esp -= 1.0;
          }
        }
      }
      else
      {
        for (int y = sy0; y <= sy1; ++y)
        {
          fill_pixel(target, x, y, color);
          esp += m;
          if (esp < -0.5)
          {
            x--;
            esp += 1.0;
          }
        }
      }
    }
  }

  void SoftwareRendererImp::rasterize_triangle(float x0, float y0,
                                               float x1, float y1,
                                               float x2, float y2,
                                               Color color)
  {
    // Task 3:
    // Implement triangle rasterization
    unsigned char* target;
    size_t width, height;
    if(sample_rate != 1) {
      x0 *= sample_rate;
      y0 *= sample_rate;
      x1 *= sample_rate;
      y1 *= sample_rate;
      x2 *= sample_rate;
      y2 *= sample_rate;
      target = sample_target;
      width = target_w * sample_rate;
      height = target_h * sample_rate;
    } else {
      target = render_target;
      width = target_w;
      height = target_h;
    }

    float max_x, max_y, middle_x, middle_y, min_x, min_y;
  
    if (y0 >= y1 && y0 >= y2)
    {
      max_x = x0;
      max_y = y0;
      if (y1 >= y2)
      {
        middle_x = x1;
        middle_y = y1;
        min_x = x2;
        min_y = y2;
      }
      else
      {
        middle_x = x2;
        middle_y = y2;
        min_x = x1;
        min_y = y1;
      }
    }
    else if (y1 >= y2)
    {
      max_x = x1;
      max_y = y1;
      if (y0 >= y2)
      {
        middle_x = x0;
        middle_y = y0;
        min_x = x2;
        min_y = y2;
      }
      else
      {
        middle_x = x2;
        middle_y = y2;
        min_x = x0;
        min_y = y0;
      }
    }
    else
    {
      max_x = x2;
      max_y = y2;
      if (y0 >= y1)
      {
        middle_x = x0;
        middle_y = y0;
        min_x = x1;
        min_y = y1;
      }
      else
      {
        middle_x = x1;
        middle_y = y1;
        min_x = x0;
        min_y = y0;
      }
    }

    float intersection_t = (middle_y - min_y) / (max_y - min_y);
    float left_x = (1.0 - intersection_t) * min_x + intersection_t * max_x;
    float left_y = (1.0 - intersection_t) * min_y + intersection_t * max_y;
    float right_x = middle_x;
    float right_y = middle_y;
    if (right_x < left_x)
    {
      swap(left_x, right_x);
      swap(left_y, right_y);
    }

    int s_min_x, s_max_x;
    // upper triangle
    float left_dxdy = (max_x - left_x) / (max_y - left_y);
    float right_dxdy = (max_x - right_x) / (max_y - right_y);
    int s_min_y = floor(middle_y);
    int s_max_y = floor(max_y);
    float temp_left_x = left_x;
    float temp_right_x = right_x;
    for (int y = s_min_y; y <= s_max_y; ++y)
    {
      s_min_x = floor(temp_left_x);
      s_max_x = floor(temp_right_x);
      for (int x = s_min_x; x <= s_max_x; ++x)
      {
        if (!(x >= 0 && x < width && y >= 0 && y < height))
        {
          continue;
        }
        fill_sample(target, x, y, color);
      }
      temp_left_x += left_dxdy;
      temp_right_x += right_dxdy;
    }

    // lower triangle
    left_dxdy = (min_x - left_x) / (min_y - left_y);
    right_dxdy = (min_x - right_x) / (min_y - right_y);
    s_min_y = floor(min_y);
    s_max_y = floor(middle_y);
    temp_left_x = left_x;
    temp_right_x = right_x;
    for (int y = s_max_y; y >= s_min_y; --y)
    {
      s_min_x = floor(temp_left_x);
      s_max_x = floor(temp_right_x);
      for (int x = s_min_x; x <= s_max_x; ++x)
      {
        if (!(x >= 0 && x < width && y >= 0 && y < height))
        {
          continue;
        }
        fill_sample(target, x, y, color);
      }
      temp_left_x -= left_dxdy;
      temp_right_x -= right_dxdy;
    }
  }

  void SoftwareRendererImp::rasterize_image(float x0, float y0,
                                            float x1, float y1,
                                            Texture &tex)
  {
    // Task 6:
    // Implement image rasterization
    unsigned char* target;
    if(sample_rate != 1) {
      target = sample_target;
      x0 *= sample_rate;
      y0 *= sample_rate;
      x1 *= sample_rate;
      y1 *= sample_rate;
    } else {
      target = render_target;
    }
    size_t min_x = fmax(0.0, x0), min_y = fmax(0.0, y0), max_x = fmin(this->target_w, x1), max_y = fmin(this->target_h, y1);
    float u, v;
    float u_scale = (x1 - x0) / tex.width;
    float v_scale = (y1 - y0) / tex.height;
    for(size_t y = min_y; y <= max_y; ++y) {
      v = (y - y0) / (y1 - y0);
      for(size_t x = min_x; x <= max_x; ++x) {
        u = (x - x0) / (x1 - x0);
        Color color = sampler->sample_trilinear(tex, u, v, u_scale, v_scale);
        fill_sample(target, x, y, color);
      }
    }
  }

  // resolve samples to render target
  void SoftwareRendererImp::resolve(void)
  {

    // Task 4:
    // Implement supersampling
    // You may also need to modify other functions marked with "Task 4".
    if(sample_rate == 1) {
      return;
    }
    size_t sample_w = target_w * sample_rate;
    size_t sample_x = 0, sample_y = 0;
    float sample_rate_square_reciprocal = 1.0 / (sample_rate * sample_rate);
    for(size_t y = 0; y < target_h; ++y) {
      sample_x = 0;
      for(size_t x = 0; x < target_w; ++x) {
        unsigned int total[4] = {0, 0, 0, 0};
        for(size_t i = 0; i < sample_rate; ++i) {
          for(size_t j = 0; j < sample_rate; ++j) {
            unsigned int p = 4 * (sample_x + j + (sample_y + i) * sample_w);
            total[0] += sample_target[p];
            total[1] += sample_target[p + 1];
            total[2] += sample_target[p + 2];
            total[3] += sample_target[p + 3];
          }
        }
        unsigned int p = 4 * (x + y * target_w);
        render_target[p] = (uint8_t) (total[0] * sample_rate_square_reciprocal);
        render_target[p + 1] = (uint8_t) (total[1] * sample_rate_square_reciprocal);
        render_target[p + 2] = (uint8_t) (total[2] * sample_rate_square_reciprocal);
        render_target[p + 3] = (uint8_t) (total[3] * sample_rate_square_reciprocal);
        sample_x += sample_rate;
      }
      sample_y += sample_rate;
    }
  }

  void SoftwareRendererImp::fill_sample(unsigned char* target, size_t sx, size_t sy, const Color& c) {
    size_t width = target_w * sample_rate;
    unsigned int p = 4 * (sx + sy * width);
    target[p] = (uint8_t)(c.r * 255 + target[p] * (1.0 - c.a));
    target[p + 1] = (uint8_t)(c.g * 255 + target[p + 1] * (1.0 - c.a));
    target[p + 2] = (uint8_t)(c.b * 255 + target[p + 2] * (1.0 - c.a));
    target[p + 3] = (uint8_t)(255 - (255 - target[p + 3]) * (1.0 - c.a));
  }

  void SoftwareRendererImp::fill_pixel(unsigned char* target, size_t x, size_t y, const Color& c) {
    size_t sample_x = x * sample_rate;
    size_t sample_y = y * sample_rate;
    size_t width = target_w * sample_rate;
    float r = c.r * 255;
    float g = c.g * 255;
    float b = c.b * 255;
    float a = c.a * 255;
    float alpha = c.a;
    for(size_t i = 0; i < sample_rate; ++i) {
      for(size_t j = 0; j < sample_rate; ++j) {
        size_t p = 4 * (sample_x + j + width * (sample_y + i));
        target[p] = uint8_t(r + target[p] * (1.0 - alpha));
        target[p + 1] = uint8_t(g + target[p + 1] * (1.0 - alpha));
        target[p + 2] = uint8_t(b + target[p + 2] * (1.0 - alpha));
        target[p + 3] = uint8_t(255 - (255 - target[p + 3]) * (1.0 - alpha));
      }
    }
  }
} // namespace CMU462

// Python - C extension for polygon module.
//
// To compile (Python2.3):
//  gcc -c polygon_ext.c -I/usr/include/python2.3 -o polygon_ext.o -Wall -O
//  gcc -shared polygon_ext.o  -o polygon_ext.so
//
// See the module polygon.py
//
//
// Ole Nielsen, GA 2004
//
// NOTE: We use anuga_int* for numeric arrays as this will work both
//       for 64 as well as 32 bit systems. See anuga_typedefs.h

#include "math.h"
#include "stdint.h"
#include "stdio.h"
#include <stdlib.h>
#include "anuga_typedefs.h"

/* Defined below; used by the segment-triangle predicate for a zero-length
 * segment. */
anuga_int __is_inside_triangle(double *point,
                             double *triangle,
                             anuga_int closed,
                             double rtol,
                             double atol);
#define YES 1
#define NO 0

inline double dist(const double x,
                   const double y)
{

  return sqrt(x * x + y * y);
}

anuga_int __point_on_line(const double x, const double y,
                        const double x0, const double y0,
                        const double x1, const double y1,
                        const double rtol,
                        const double atol)
{
  /*Determine whether a point is on a line segment

    Input: x, y, x0, x0, x1, y1: where
        point is given by x, y
  line is given by (x0, y0) and (x1, y1)

  */

  double a0, a1, a_normal0, a_normal1, b0, b1, len_a, len_b;
  double a_dot_b;
  double nominator, denominator;
  anuga_int is_parallel;

  a0 = x - x0;
  a1 = y - y0;

  a_normal0 = a1;
  a_normal1 = -a0;

  b0 = x1 - x0;
  b1 = y1 - y0;

  nominator = fabs(a_normal0 * b0 + a_normal1 * b1);
  denominator = b0 * b0 + b1 * b1;

  // Determine if line is parallel to point vector up to a tolerance
  is_parallel = 0;
  if (denominator == 0.0)
  {
    // Use absolute tolerance
    if (nominator <= atol)
    {
      is_parallel = 1;
    }
  }
  else
  {
    // Denominator is positive - use relative tolerance
    if (nominator / denominator <= rtol)
    {
      is_parallel = 1;
    }
  }

  if (is_parallel)
  {
    // Point is somewhere on the infinite extension of the line
    // subject to specified absolute tolerance

    len_a = dist(a0, a1); // sqrt(a0*a0 + a1*a1);
    len_b = dist(b0, b1); // sqrt(b0*b0 + b1*b1);

    a_dot_b = a0 * b0 + a1 * b1;

    if (a_dot_b >= -1.0e-308 && -1.0e-15 <= len_b - len_a)
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }
  else
  {
    return 0;
  }
}

//  public domain function by Darel Rex Finley, 2006
//  http://www.alienryderflex.com/intersect/
//
//  Determines the intersection point of the line segment defined by points A and B
//  with the line segment defined by points C and D.
//
//  Returns YES if the intersection point was found, and stores that point in X,Y.
//  Returns NO if there is no determinable intersection point, in which case X,Y will
//  be unmodified.

anuga_int __line_segment_intersection(
    double Ax, double Ay,
    double Bx, double By,
    double Cx, double Cy,
    double Dx, double Dy,
    double *X, double *Y)
{

  double distAB, theCos, theSin, newX, ABpos;

  //  Fail if either line segment is zero-length.
  if ((Ax == Bx && Ay == By) || (Cx == Dx && Cy == Dy))
    return NO;

  //  Fail if the segments share an end-point.
  if ((Ax == Cx && Ay == Cy) || (Bx == Cx && By == Cy) || (Ax == Dx && Ay == Dy) || (Bx == Dx && By == Dy))
  {
    return NO;
  }

  //  (1) Translate the system so that point A is on the origin.
  Bx -= Ax;
  By -= Ay;
  Cx -= Ax;
  Cy -= Ay;
  Dx -= Ax;
  Dy -= Ay;

  //  Discover the length of segment A-B.
  distAB = sqrt(Bx * Bx + By * By);

  //  (2) Rotate the system so that point B is on the positive X axis.
  theCos = Bx / distAB;
  theSin = By / distAB;
  newX = Cx * theCos + Cy * theSin;
  Cy = Cy * theCos - Cx * theSin;
  Cx = newX;
  newX = Dx * theCos + Dy * theSin;
  Dy = Dy * theCos - Dx * theSin;
  Dx = newX;

  //  Fail if segment C-D doesn't cross line A-B.
  if ((Cy < 0. && Dy < 0.) || (Cy >= 0. && Dy >= 0.))
    return NO;

  //  (3) Discover the position of the intersection point along line A-B.
  ABpos = Dx + (Cx - Dx) * Dy / (Dy - Cy);

  //  Fail if segment C-D crosses line A-B outside of segment A-B.
  if (ABpos < 0. || ABpos > distAB)
    return NO;

  //  (4) Apply the discovered position to line A-B in the original coordinate system.
  *X = Ax + ABpos * theCos;
  *Y = Ay + ABpos * theSin;

  //  Success.
  return YES;
}

/*
WORK IN PROGRESS TO OPTIMISE INTERSECTION
anuga_int __intersection(double x0, double y0,
       double x1, double y1) {


    x0 = line0[0,0]; y0 = line0[0,1]
    x1 = line0[1,0]; y1 = line0[1,1]

    x2 = line1[0,0]; y2 = line1[0,1]
    x3 = line1[1,0]; y3 = line1[1,1]

    denom = (y3-y2)*(x1-x0) - (x3-x2)*(y1-y0)
    u0 = (x3-x2)*(y0-y2) - (y3-y2)*(x0-x2)
    u1 = (x2-x0)*(y1-y0) - (y2-y0)*(x1-x0)

    if allclose(denom, 0.0):
        # Lines are parallel - check if they coincide on a shared a segment

        if allclose( [u0, u1], 0.0 ):
            # We now know that the lines if continued coincide
            # The remaining check will establish if the finite lines share a segment

            line0_starts_on_line1 = line0_ends_on_line1 =\
            line1_starts_on_line0 = line1_ends_on_line0 = False

            if point_on_line([x0, y0], line1):
                line0_starts_on_line1 = True

            if point_on_line([x1, y1], line1):
                line0_ends_on_line1 = True

            if point_on_line([x2, y2], line0):
                line1_starts_on_line0 = True

            if point_on_line([x3, y3], line0):
                line1_ends_on_line0 = True

            if not(line0_starts_on_line1 or line0_ends_on_line1\
               or line1_starts_on_line0 or line1_ends_on_line0):
                # Lines are parallel and would coincide if extended, but not as they are.
                return 3, None


            # One line fully included in the other. Use direction of included line
            if line0_starts_on_line1 and line0_ends_on_line1:
                # Shared segment is line0 fully included in line1
                segment = array([[x0, y0], [x1, y1]])

            if line1_starts_on_line0 and line1_ends_on_line0:
                # Shared segment is line1 fully included in line0
                segment = array([[x2, y2], [x3, y3]])


            # Overlap with lines are oriented the same way
            if line0_starts_on_line1 and line1_ends_on_line0:
                # Shared segment from line0 start to line 1 end
                segment = array([[x0, y0], [x3, y3]])

            if line1_starts_on_line0 and line0_ends_on_line1:
                # Shared segment from line1 start to line 0 end
                segment = array([[x2, y2], [x1, y1]])


            # Overlap in opposite directions - use direction of line0
            if line0_starts_on_line1 and line1_starts_on_line0:
                # Shared segment from line0 start to line 1 end
                segment = array([[x0, y0], [x2, y2]])

            if line0_ends_on_line1 and line1_ends_on_line0:
                # Shared segment from line0 start to line 1 end
                segment = array([[x3, y3], [x1, y1]])


            return 2, segment
        else:
            # Lines are parallel but they do not coincide
            return 4, None #FIXME (Ole): Add distance here instead of None

    else:
        # Lines are not parallel or coinciding
        u0 = u0/denom
        u1 = u1/denom

        x = x0 + u0*(x1-x0)
        y = y0 + u0*(y1-y0)

        # Sanity check - can be removed to speed up if needed
        assert allclose(x, x2 + u1*(x3-x2))
        assert allclose(y, y2 + u1*(y3-y2))

        # Check if point found lies within given line segments
        if 0.0 <= u0 <= 1.0 and 0.0 <= u1 <= 1.0:
            # We have intersection

            return 1, array([x, y])
        else:
            # No intersection
            return 0, None


}
*/

anuga_int __interpolate_polyline(anuga_int number_of_nodes,
                               anuga_int number_of_points,
                               double *data,
                               double *polyline_nodes,
                               anuga_int *gauge_neighbour_id,
                               double *interpolation_points,
                               double *interpolated_values,
                               double rtol,
                               double atol)
{

  anuga_int j, i, neighbour_id;
  double x0, y0, x1, y1, x, y;
  double segment_len, segment_delta, slope, alpha;

  for (j = 0; j < number_of_nodes; j++)
  {

    neighbour_id = gauge_neighbour_id[j];

    // gauge_neighbour_id[j] names the node that segment j runs to, and is -1
    // for the last node of a polyline. It is passed in rather than assumed to
    // be j+1 so a set of disjoint polylines can share one array; keep it.

    if (neighbour_id >= 0)
    {
      x0 = polyline_nodes[2 * j];
      y0 = polyline_nodes[2 * j + 1];

      x1 = polyline_nodes[2 * neighbour_id];
      y1 = polyline_nodes[2 * neighbour_id + 1];

      segment_len = dist(x1 - x0, y1 - y0);
      segment_delta = data[neighbour_id] - data[j];
      slope = segment_delta / segment_len;

      for (i = 0; i < number_of_points; i++)
      {
        x = interpolation_points[2 * i];
        y = interpolation_points[2 * i + 1];

        if (__point_on_line(x, y, x0, y0, x1, y1, rtol, atol))
        {
          alpha = dist(x - x0, y - y0);
          interpolated_values[i] = slope * alpha + data[j];
        }
      }
    }
  }

  return 0;
}

anuga_int __triangle_polygon_overlap(double *polygon,
                                   double *triangle,
                                   anuga_int polygon_number_of_vertices)
{
  anuga_int i, ii, j, jj, A, B;
  double p0_x, p0_y, p1_x, p1_y, pp_x, pp_y;
  double t0_x, t0_y, t1_x, t1_y, tp_x, tp_y;
  double u_x, u_y, v_x, v_y, w_x, w_y;
  double u_dot_tp, v_dot_tp, v_dot_pp, w_dot_pp;
  double a, b;

  p0_x = polygon[0];
  p0_y = polygon[1];

  A = 0;
  B = 0;

  for (i = 1; i < polygon_number_of_vertices + 1; i++)
  {
    ii = i % polygon_number_of_vertices;

    p1_x = polygon[2 * ii];
    p1_y = polygon[2 * ii + 1];

    pp_x = -(p1_y - p0_y);
    pp_y = p1_x - p0_x;

    t0_x = triangle[0];
    t0_y = triangle[1];

    for (j = 1; j < 4; j++)
    {
      jj = j % 3;

      t1_x = triangle[2 * jj];
      t1_y = triangle[2 * jj + 1];

      tp_x = -(t1_y - t0_y); // perpendicular to triangle vector
      tp_y = t1_x - t0_x;    // perpendicular to polygon vector

      u_x = p1_x - p0_x;
      u_y = p1_y - p0_y;
      v_x = t0_x - p0_x;
      v_y = t0_y - p0_y;
      w_x = t1_x - t0_x;
      w_y = t1_y - t0_y;

      u_dot_tp = (u_x * tp_x) + (u_y * tp_y);

      if (u_dot_tp != 0.0f) // Vectors are not parallel
      {
        v_dot_tp = (v_x * tp_x) + (v_y * tp_y);
        v_dot_pp = (v_x * pp_x) + (v_y * pp_y);
        w_dot_pp = (w_x * pp_x) + (w_y * pp_y);

        a = v_dot_tp / u_dot_tp;
        b = -v_dot_pp / w_dot_pp;

        if (a >= 0.0f && a <= 1.0f && b >= 0.0f && b <= 1.0f)
        {
          return 1; // overlap
        }

        if (b >= 0.0f && b <= 1.0f && a > 1.0f)
        {
          A++;
        }

        if (a >= 0.0f && a <= 1.0f && b > 1.0f)
        {
          B++;
        }

        if (A == 4 || B == 3)
        {
          return 1; // overlap
        }
      }

      t0_x = t1_x;
      t0_y = t1_y;
    }

    p0_x = p1_x;
    p0_y = p1_y;
  }

  return 0; // no overlap
}

anuga_int __polygon_overlap(double *polygon,
                          double *triangles,
                          anuga_int *indices,
                          anuga_int M, // number of triangles
                          anuga_int polygon_number_of_vertices)
{
  double *triangle;
  anuga_int i, inside_index, outside_index;

  inside_index = 0;      // Keep track of triangles that overlap
  outside_index = M - 1; // Keep track of triangles that don't overlap (starting from end)

  for (i = 0; i < M; i++)
  {
    triangle = triangles + 6 * i;

    if (__triangle_polygon_overlap(polygon,
                                   triangle,
                                   polygon_number_of_vertices))
    {
      indices[inside_index] = i;
      inside_index++;
    }
    else
    {
      indices[outside_index] = i;
      outside_index -= 1;
    }
  }

  return inside_index;
}

/* Does the segment p0->p1 pass THROUGH or ALONG the closed triangle?
 *
 * The segment is clipped against the triangle's three half-planes
 * (Cyrus-Beck) to the parameter interval [amin, amax] of the part of the
 * segment inside the triangle, and the triangle counts as intersected when
 * that part has positive length. So:
 *
 *   - a segment crossing the interior, ending inside, or lying along an
 *     edge counts (a segment along a shared edge counts for the triangles
 *     on BOTH sides -- the selection is symmetric, a property of the
 *     geometry);
 *   - a segment that only touches a vertex, or ends exactly on an edge from
 *     outside, does not: the contact has zero length.
 *
 * Tolerances are lengths relative to the size of the segment and triangle,
 * so a segment within rounding of a mesh line is treated as on it whichever
 * way the last bit fell. The previous test compared exact intersection
 * parameters against 0 and 1 and counted point contacts, so triangles
 * touched at a vertex or along an edge were selected or not according to
 * floating-point contraction and evaluation order -- the same script gave
 * a 68% larger inlet area under one compiler than another when its exchange
 * line lay on mesh lines (anuga-community/anuga_core#231). */
anuga_int __triangle_line_intersect(double *line,
                                  double *triangle)
{
  const double p0_x = line[0], p0_y = line[1];
  const double p1_x = line[2], p1_y = line[3];
  const double u_x = p1_x - p0_x, u_y = p1_y - p0_y;
  const double u_len = sqrt(u_x * u_x + u_y * u_y);

  /* A length scale for the tolerances: the segment and the triangle. */
  double scale = u_len;
  double area2 = 0.0;
  anuga_int i;
  for (i = 0; i < 3; i++)
  {
    const anuga_int k = (i + 1) % 3;
    const double e_x = triangle[2 * k] - triangle[2 * i];
    const double e_y = triangle[2 * k + 1] - triangle[2 * i + 1];
    const double e_len = sqrt(e_x * e_x + e_y * e_y);
    if (e_len > scale) scale = e_len;
  }
  const double tol = 1.0e-10 * scale;           /* a length */

  area2 = (triangle[2] - triangle[0]) * (triangle[5] - triangle[1])
        - (triangle[4] - triangle[0]) * (triangle[3] - triangle[1]);
  if (fabs(area2) <= tol * tol)
  {
    return 0;                                   /* degenerate triangle */
  }
  const double orient = (area2 > 0.0) ? 1.0 : -1.0;

  if (u_len <= tol)
  {
    /* A zero-length segment has no interior and selects nothing. It is
     * what the two side edges of a structure's apron polygon are when
     * apron = 0, and treating it as a point test pulled in the triangle
     * across the channel bank at each end of a bridge's exchange line,
     * which moved the bridge's discharge (the HEC-RAS bridge regressions
     * caught it). A caller that wants a point test has is_inside_triangle. */
    return 0;
  }

  double amin = 0.0, amax = 1.0;
  for (i = 0; i < 3; i++)
  {
    const anuga_int k = (i + 1) % 3;
    const double ta_x = triangle[2 * i], ta_y = triangle[2 * i + 1];
    const double e_x = triangle[2 * k] - ta_x;
    const double e_y = triangle[2 * k + 1] - ta_y;
    const double e_len = sqrt(e_x * e_x + e_y * e_y);
    if (e_len <= tol)
    {
      continue;                                 /* zero-length edge */
    }
    /* Unit inward normal of this edge (left of the edge for a
     * counter-clockwise triangle), so d0 and du are lengths. */
    const double n_x = orient * (-e_y) / e_len;
    const double n_y = orient * e_x / e_len;
    const double d0 = n_x * (p0_x - ta_x) + n_y * (p0_y - ta_y);
    const double du = n_x * u_x + n_y * u_y;

    if (fabs(du) <= tol)
    {
      /* Segment parallel to this edge: it is inside this half-plane
       * everywhere or nowhere. Within tolerance of the edge counts as on
       * it, on either side. */
      if (d0 < -tol)
      {
        return 0;
      }
      continue;
    }
    const double a = -d0 / du;
    if (du > 0.0)
    {
      if (a > amin) amin = a;                   /* entering */
    }
    else
    {
      if (a < amax) amax = a;                   /* leaving */
    }
  }

  /* Positive length inside the triangle, not a point contact. */
  return ((amax - amin) * u_len > tol) ? 1 : 0;
}

anuga_int __line_intersect(double *line,
                         double *triangles,
                         anuga_int *indices,
                         anuga_int M) // number of triangles
{
  double *triangle;
  anuga_int i, inside_index, outside_index;

  inside_index = 0;      // Keep track of triangles that intersect
  outside_index = M - 1; // Keep track of triangles that don't intersect (starting from end)

  for (i = 0; i < M; i++)
  {
    triangle = triangles + 6 * i;

    if (__triangle_line_intersect(line,
                                  triangle))
    {
      indices[inside_index] = i;
      inside_index++;
    }
    else
    {
      indices[outside_index] = i;
      outside_index -= 1;
    }
  }

  return inside_index;
}

anuga_int __is_inside_triangle(double *point,
                             double *triangle,
                             anuga_int closed,
                             double rtol,
                             double atol)
{

  double vx, vy, v0x, v0y, v1x, v1y;
  double a00, a10, a01, a11, b0, b1;
  double denom, alpha, beta;

  double x, y; // Point coordinates
  anuga_int i, j, res;

  x = point[0];
  y = point[1];

  // Quickly reject points that are clearly outside
  if ((x < triangle[0]) &&
      (x < triangle[2]) &&
      (x < triangle[4]))
    return 0;

  if ((x > triangle[0]) &&
      (x > triangle[2]) &&
      (x > triangle[4]))
    return 0;

  if ((y < triangle[1]) &&
      (y < triangle[3]) &&
      (y < triangle[5]))
    return 0;

  if ((y > triangle[1]) &&
      (y > triangle[3]) &&
      (y > triangle[5]))
    return 0;

  // v0 = C-A
  v0x = triangle[4] - triangle[0];
  v0y = triangle[5] - triangle[1];

  // v1 = B-A
  v1x = triangle[2] - triangle[0];
  v1y = triangle[3] - triangle[1];

  // First check if point lies wholly inside triangle
  a00 = v0x * v0x + v0y * v0y; // innerproduct(v0, v0)
  a01 = v0x * v1x + v0y * v1y; // innerproduct(v0, v1)
  a10 = a01;                   // innerproduct(v1, v0)
  a11 = v1x * v1x + v1y * v1y; // innerproduct(v1, v1)

  denom = a11 * a00 - a01 * a10;

  if (fabs(denom) > 0.0)
  {
    // v = point-A
    vx = x - triangle[0];
    vy = y - triangle[1];

    b0 = v0x * vx + v0y * vy; // innerproduct(v0, v)
    b1 = v1x * vx + v1y * vy; // innerproduct(v1, v)

    alpha = (b0 * a11 - b1 * a01) / denom;
    beta = (b1 * a00 - b0 * a10) / denom;

    if ((alpha > 0.0) && (beta > 0.0) && (alpha + beta < 1.0))
      return 1;
  }

  if (closed)
  {
    // Check if point lies on one of the edges

    for (i = 0; i < 3; i++)
    {
      j = (i + 1) % 3; // Circular index into triangle array
      res = __point_on_line(x, y,
                            triangle[2 * i], triangle[2 * i + 1],
                            triangle[2 * j], triangle[2 * j + 1],
                            rtol, atol);
      if (res)
        return 1;
    }
  }

  // Default return if point is outside triangle
  return 0;
}

anuga_int __separate_points_by_polygon(const anuga_int M, // Number of points
                                     const anuga_int N, // Number of polygon vertices
                                     double *points,
                                     double *polygon,
                                     anuga_int *indices, // M-Array for storage indices
                                     const anuga_int closed,
                                     const anuga_int verbose,
                                     const double rtol,  // tolerances for the
                                     const double atol)  // point-on-edge test
{

  double minpx, maxpx, minpy, maxpy;
  anuga_int outside_index, inside_index;

  // Find min and max of poly used for optimisation when points
  // are far away from polygon

  minpx = polygon[0];
  maxpx = minpx;
  minpy = polygon[1];
  maxpy = minpy;
#pragma omp parallel for reduction(min : minpx, minpy) reduction(max : maxpx, maxpy)
  for (int i = 0; i < N; i++)
  {
    double px_i = polygon[2 * i];
    double py_i = polygon[2 * i + 1];

    if (px_i < minpx)
      minpx = px_i;
    if (px_i > maxpx)
      maxpx = px_i;
    if (py_i < minpy)
      minpy = py_i;
    if (py_i > maxpy)
      maxpy = py_i;
  }

  // Main loop, in two passes. The inside/outside test for each point is
  // independent, so it runs in parallel into a flag array; the compaction
  // into `indices` is order dependent (inside points in increasing k from
  // the front, outside points in decreasing k from the back) and stays
  // serial, so the output is identical to the old single loop.
  unsigned char *inside_flag = (unsigned char *)malloc((size_t)(M > 0 ? M : 1));
  if (inside_flag == NULL)
    return -1;

#pragma omp parallel for schedule(static)
  for (int k = 0; k < M; k++)
  {
    double x = points[2 * k];
    double y = points[2 * k + 1];

    int inside = 0;

    // Optimisation
    if ((x > maxpx) || (x < minpx) || (y > maxpy) || (y < minpy))
    {
      // Nothing
    }
    else
    {
      // Check polygon
      for (int i = 0; i < N; i++)
      {
        int j = (i + 1) % N;

        double px_i = polygon[2 * i];
        double py_i = polygon[2 * i + 1];
        double px_j = polygon[2 * j];
        double py_j = polygon[2 * j + 1];

        // Check for case where point is contained in line segment
        if (__point_on_line(x, y, px_i, py_i, px_j, py_j, rtol, atol))
        {
          if (closed == 1)
          {
            inside = 1;
          }
          else
          {
            inside = 0;
          }
          break;
        }
        else
        {
          // Check if truly inside polygon
          if (((py_i < y) && (py_j >= y)) ||
              ((py_j < y) && (py_i >= y)))
          {
            if (px_i + (y - py_i) / (py_j - py_i) * (px_j - px_i) < x)
              inside = 1 - inside;
          }
        }
      }
    }
    inside_flag[k] = (unsigned char)inside;
  } // End k (parallel)

  inside_index = 0;      // Keep track of points inside
  outside_index = M - 1; // Keep track of points outside (starting from end)
  for (int k = 0; k < M; k++)
  {
    if (inside_flag[k])
    {
      indices[inside_index] = k;
      inside_index += 1;
    }
    else
    {
      indices[outside_index] = k;
      outside_index -= 1;
    }
  }
  free(inside_flag);

  return inside_index;
}

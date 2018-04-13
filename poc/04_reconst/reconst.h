#include "vect.hpp"

using Scal = double;
using Vect = geom::GVect<Scal, 3>;


// line:
// a: volume fraction 
// s: area of fluid 

/**
 * gfs_line_s:
 * @m: normal to the line.
 * @a: line constant.
 *
 * Returns: the area of the fraction of s cell lying under the line
 * (@m,@a).
 */
Scal gfs_line_s (const Vect& m, Scal a);

/**
 * gfs_line_a:
 * @m: s #Vect.
 * @c: s volume fraction.
 *
 * Returns: the value @a such that the area of s square cell
 * lying under the line defined by @m.@x = @a is equal to @c. 
 */
Scal gfs_line_a (const Vect& m, Scal c);

/**
 * gfs_line_c:
 * @m: normal to the line.
 * @a: line constant.
 * @s: area of cell fraction.
 * @p: s #Vect.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * s square cell lying under the line (@m,@a).
 */
void gfs_line_c (const Vect& m, Scal a, Scal s, Vect& p);

/**
 * gfs_line_sc:
 * @m: normal to the line.
 * @a: line constant.
 * @p: s #Vect.
 *
 * Fills @p with the position of the center of area of the fraction of
 * s square cell lying under the line (@m,@a).
 *
 * Returns: the length of the facet.
 */
Scal gfs_line_sc (const Vect& m, Scal a, Vect& p);

/**
 * gfs_plane_volume:
 * @m: normal to the plane.
 * @a: plane constant.
 *
 * Returns: the volume of s cell lying under the plane (@m,@a).
 */
Scal gfs_plane_volume (const Vect& m, Scal a);

/**
 * gfs_plane_a:
 * @m: s #Vect.
 * @c: s volume fraction.
 *
 * Returns: the value @a such that the volume of s cubic cell
 * lying under the plane defined by @m.@x = @a is equal to @c. 
 */
Scal gfs_plane_a (const Vect& m, Scal c);

/**
 * gfs_planec:
 * @m: normal to the plane.
 * @a: plane constant.
 * @s: volume of cell fraction.
 * @p: s #Vect.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * s cubic cell lying under the plane (@m,@a).
 */
void gfs_planec (const Vect& m, Scal a, Scal s, Vect& p);

/**
 * gfs_plane_sc:
 * @m: normal to the plane.
 * @a: plane constant.
 * @p: s #Vect.
 *
 * Fills @p with the position of the center of area of the fraction of
 * s cubic cell lying under the plane (@m,@a).
 *
 * Returns: the area of the facet.
 */
Scal gfs_plane_sc (const Vect& m, Scal a, Vect& p);


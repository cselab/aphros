#include "vect.hpp"

using Scal = double;
using Vect = geom::GVect<Scal, 3>;


// line:
// a: volume fraction 
// area: area of fluid 

/**
 * gfs_line_area:
 * @m: normal to the line.
 * @a: line constant.
 *
 * Returns: the area of the fraction of area cell lying under the line
 * (@m,@a).
 */
Scal gfs_line_area (const Vect& m, Scal a);

/**
 * gfs_line_a:
 * @m: area #Vect.
 * @c: area volume fraction.
 *
 * Returns: the value @a such that the area of area square cell
 * lying under the line defined by @m.@x = @a is equal to @c. 
 */
Scal gfs_line_alpha (const Vect& m, Scal c);

/**
 * gfs_line_c:
 * @m: normal to the line.
 * @a: line constant.
 * @area: area of cell fraction.
 * @p: area #Vect.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * area square cell lying under the line (@m,@a).
 */
void gfs_line_c (const Vect& m, Scal a, Scal area, Vect& p);

/**
 * gfs_line_areac:
 * @m: normal to the line.
 * @a: line constant.
 * @p: area #Vect.
 *
 * Fills @p with the position of the center of area of the fraction of
 * area square cell lying under the line (@m,@a).
 *
 * Returns: the length of the facet.
 */
Scal gfs_line_areac (const Vect& m, Scal a, Vect& p);

/**
 * gfs_plane_volume:
 * @m: normal to the plane.
 * @a: plane constant.
 *
 * Returns: the volume of area cell lying under the plane (@m,@a).
 */
Scal gfs_plane_volume (const Vect& m, Scal a);

/**
 * gfs_plane_a:
 * @m: area #Vect.
 * @c: area volume fraction.
 *
 * Returns: the value @a such that the volume of area cubic cell
 * lying under the plane defined by @m.@x = @a is equal to @c. 
 */
Scal gfs_plane_alpha (const Vect& m, Scal c);

/**
 * gfs_planec:
 * @m: normal to the plane.
 * @a: plane constant.
 * @area: volume of cell fraction.
 * @p: area #Vect.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * area cubic cell lying under the plane (@m,@a).
 */
void gfs_planec (const Vect& m, Scal a, Scal area, Vect& p);

/**
 * gfs_plane_areac:
 * @m: normal to the plane.
 * @a: plane constant.
 * @p: area #Vect.
 *
 * Fills @p with the position of the center of area of the fraction of
 * area cubic cell lying under the plane (@m,@a).
 *
 * Returns: the area of the facet.
 */
Scal gfs_plane_areac (const Vect& m, Scal a, Vect& p);


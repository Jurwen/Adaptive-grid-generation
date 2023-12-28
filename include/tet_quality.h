//
//  tet_quality.h
//  adaptive_mesh_refinement
//
//  Created by Yiwen Ju on 12/28/23.
//  the measure of radius ratio is following this git: https://github.com/sandialabs/verdict/tree/master

#ifndef tet_quality_h
#define tet_quality_h

double tet_radius_ratio(int /*num_nodes*/, const double coordinates[][3])
{

  // Determine side vectors
  VerdictVector side[6];

  side[0].set(coordinates[1][0] - coordinates[0][0], coordinates[1][1] - coordinates[0][1],
    coordinates[1][2] - coordinates[0][2]);

  side[1].set(coordinates[2][0] - coordinates[1][0], coordinates[2][1] - coordinates[1][1],
    coordinates[2][2] - coordinates[1][2]);

  side[2].set(coordinates[0][0] - coordinates[2][0], coordinates[0][1] - coordinates[2][1],
    coordinates[0][2] - coordinates[2][2]);

  side[3].set(coordinates[3][0] - coordinates[0][0], coordinates[3][1] - coordinates[0][1],
    coordinates[3][2] - coordinates[0][2]);

  side[4].set(coordinates[3][0] - coordinates[1][0], coordinates[3][1] - coordinates[1][1],
    coordinates[3][2] - coordinates[1][2]);

  side[5].set(coordinates[3][0] - coordinates[2][0], coordinates[3][1] - coordinates[2][1],
    coordinates[3][2] - coordinates[2][2]);

  VerdictVector numerator = side[3].length_squared() * (side[2] * side[0]) +
    side[2].length_squared() * (side[3] * side[0]) + side[0].length_squared() * (side[3] * side[2]);

  double area_sum;
  area_sum = ((side[2] * side[0]).length() + (side[3] * side[0]).length() +
               (side[4] * side[1]).length() + (side[3] * side[2]).length()) *
    0.5;

  double volume = tet_volume(4, coordinates);

  if (std::abs(volume) < VERDICT_DBL_MIN)
  {
    return (double)VERDICT_DBL_MAX;
  }
  else
  {
    const double radius_ratio = numerator.length() * area_sum / (108 * volume * volume);
    return fix_range(radius_ratio);
  }
}

#endif /* tet_quality_h */

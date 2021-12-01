//
// Created by Baciu Ioana on 22.11.2021.
//

#include "imageFormationUtils.h"
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv)
{
    FILE *le_file = fopen(argv[2], "w");
    struct point3d *points;
    int N_v = 0;
    points = readOff(argv[1], &N_v);
//    printf("%d", N_v);
    centerThePCL(points, N_v);
    struct point3d *transformed_points = rigid_transformation(points, N_v, -1.57,3.2,-0.2);
    int N2 = N_v;
//    struct point3d *filtered_points = filter_points(transformed_points, &N2);
//    struct point2d *other_points = pinhole_projection(filtered_points, N2, (float )100000);
    struct point2d *other_points = another_filter_points(transformed_points,  N2, (float )1);
    image_structure_rgb_t *img = pixel_coordinate_projection(other_points, N2, 500, 500, 200,200,200, 200);

    write_image_to_file_rgb(img, le_file);
    //
    return 0;
}
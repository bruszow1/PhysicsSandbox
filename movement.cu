#include <algorithm>
#include <cmath>
#include "cublas_v2.h"
#include "nppcore.h"
#include "nppi.h"
#include <thrust/transform.h>
#include <thrust/extrema.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/execution_policy.h>
#include "movement.h"
#include "shape.h"
#include "utils.h"
#include "draw.h"

cublasStatus_t rotate_vector(double *target, double *output, double x_rotate, double y_rotate, double z_rotate, cublasHandle_t *handle) {

    const double rot_matrix[9] = { cos(x_rotate) * cos(y_rotate), sin(y_rotate), -sin(x_rotate) * cos(y_rotate),
        -cos(x_rotate) * sin(y_rotate) * cos(z_rotate) + sin(x_rotate) * sin(z_rotate), cos(y_rotate) * cos(z_rotate),
        sin(x_rotate) * sin(y_rotate) * cos(z_rotate) + cos(x_rotate) * sin(z_rotate),
        cos(x_rotate) * sin(y_rotate) * sin(z_rotate) + sin(x_rotate) * cos(z_rotate), -cos(y_rotate) * sin(z_rotate),
        -sin(x_rotate) * sin(y_rotate) * cos(z_rotate) * sin(z_rotate) + cos(x_rotate) * cos(z_rotate)};

    // alpha = 1 so matrix is unchanged, beta = 0 b/c the + C portion isn't needed
    const double alpha = 1.0;
    const double beta = 0;

    // allocate memory
    double *gpu_rot_matrix;
    cudaMalloc((void**) &gpu_rot_matrix, 9 * sizeof(double));
    double *gpu_x;
    cudaMalloc(&gpu_x, 3 * sizeof(double));
    double *gpu_y;
    cudaMalloc(&gpu_y, 3 * sizeof(double));

    // copy matrix and vectors to GPU memory
    cudaMemcpy(gpu_rot_matrix, rot_matrix, 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_x, target, 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_y, output, 3 * sizeof(double), cudaMemcpyHostToDevice);

    cublasStatus_t status = cublasDgemv_v2(*handle, CUBLAS_OP_N, 3, 3, &alpha, gpu_rot_matrix, 3, gpu_x, 1, &beta, gpu_y, 1);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("Cublas error: %d\n", status);
    }

    // copy result from GPU memory
    cudaMemcpy(output, gpu_y, 3 * sizeof(double), cudaMemcpyDeviceToHost);

    // clean up memory
    cudaFree(gpu_rot_matrix);
    cudaFree(gpu_x);
    cudaFree(gpu_y);
    return status;
}

void shift_vector(double *origin_vector, double x_move, double y_move, double z_move)
{
    origin_vector[0] += x_move;
    origin_vector[1] += y_move;
    origin_vector[2] += z_move;
}

void rotate_vector(double target[3], double output[3], double x_rotate, double y_rotate, double z_rotate)
{
    double init_x = target[0];
    double init_y = target[1];
    double init_z = target[2];

    // Use 3d rotation matrix
    output[0] = init_x * cos(x_rotate) * cos(y_rotate) + init_y * (sin(x_rotate) * sin(z_rotate) - cos(x_rotate) * sin(y_rotate) * cos(z_rotate))
                + init_z * (sin(x_rotate) * cos(z_rotate) + cos(x_rotate) * sin(y_rotate) * sin(z_rotate));
    output[1] = init_x * sin(y_rotate) + init_y * cos(y_rotate) * cos(z_rotate) - init_z * cos(y_rotate) * sin(z_rotate);
    output[2] = -1 * init_x * sin(x_rotate) * cos(y_rotate) + init_y * (sin(x_rotate) * sin(y_rotate) * cos(z_rotate) + cos(x_rotate) * sin(z_rotate))
                + init_z * (cos(x_rotate) * cos(z_rotate) - sin(x_rotate) * sin(y_rotate) * sin(z_rotate));
}

void rotate_around(Shape *shape, double x_rotate, double y_rotate, double z_rotate, double x_fixed, double y_fixed, double z_fixed, cublasHandle_t *handle)
{
    // get baseline orientation of input coords and shape origin
    double shape_axis_vector[3] = { x_fixed - shape->origin[0], y_fixed - shape->origin[1], z_fixed - shape->origin[2] };

    const double rot_matrix[9] = { cos(x_rotate) * cos(y_rotate), sin(y_rotate), -sin(x_rotate) * cos(y_rotate),
        -cos(x_rotate) * sin(y_rotate) * cos(z_rotate) + sin(x_rotate) * sin(z_rotate), cos(y_rotate) * cos(z_rotate),
        sin(x_rotate) * sin(y_rotate) * cos(z_rotate) + cos(x_rotate) * sin(z_rotate),
        cos(x_rotate) * sin(y_rotate) * sin(z_rotate) + sin(x_rotate) * cos(z_rotate), -cos(y_rotate) * sin(z_rotate),
        -sin(x_rotate) * sin(y_rotate) * cos(z_rotate) * sin(z_rotate) + cos(x_rotate) * cos(z_rotate)};

    // alpha = 1 so matrix is unchanged, beta = 0 b/c the + C portion isn't needed
    const double alpha = 1.0;
    const double beta = 0;

    // allocate GPU memory
    unsigned int vector_count = 3 * shape->triangle_count + shape->vertex_count + 1;
    unsigned int double_count = 3 * vector_count;
    double *gpu_rot_matrix;
    cudaMalloc((void**) &gpu_rot_matrix, 9 * sizeof(double));
    double *gpu_B;
    cudaMalloc(&gpu_B, double_count * sizeof(double));
    double *gpu_C;
    cudaMalloc(&gpu_C, double_count * sizeof(double));

    // copy vectors and rotation matrix to GPU memory
    cudaMemcpy(gpu_rot_matrix, rot_matrix, 9 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_B, shape_axis_vector, 3 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_B + 3, shape->vector_array[0], (double_count - 3) * sizeof(double), cudaMemcpyHostToDevice);

    cublasStatus_t status = cublasDgemm_v2(*handle, CUBLAS_OP_N, CUBLAS_OP_N, 3, vector_count, 3,  &alpha, gpu_rot_matrix, 3, gpu_B, 3, &beta, gpu_C, 3);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("Cublas error: %d\n", status);
    }

    // copy results from GPU memory
    cudaMemcpy(shape_axis_vector, gpu_C, 3 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&shape->vector_array[0], gpu_C + 3, (double_count - 3) * sizeof(double), cudaMemcpyDeviceToHost);

    // clean up GPU memory
    cudaFree(gpu_rot_matrix);
    cudaFree(gpu_B);
    cudaFree(gpu_C);

    // shift shape so target point is at original position
    double x_correction = x_fixed - (shape_axis_vector[0] + shape->origin[0]);
    double y_correction = y_fixed - (shape_axis_vector[1] + shape->origin[1]);
    double z_correction = z_fixed - (shape_axis_vector[2] + shape->origin[2]);
    shift_vector(shape->origin, x_correction, y_correction, z_correction);
}

void adjust_movement(Shape *shape, double initial_timestep, double current_timestep)
{
    double delta_t = current_timestep - initial_timestep;
    double x_displacement = shape->velocity[0] * delta_t;
    double y_displacement = shape->velocity[1] * delta_t;
    double z_displacement = shape->velocity[2] * delta_t;

    // move shape to correct coordinates
    shift_vector(shape->origin, x_displacement, y_displacement, z_displacement);
}

void adjust_rotation(Shape *shape, double initial_timestep, double current_timestep, cublasHandle_t *handle)
{
    double delta_t = current_timestep - initial_timestep;

    double x_angular_displacement = shape->angular_velocity[0] * delta_t;
    double y_angular_displacement = shape->angular_velocity[1] * delta_t;
    double z_angular_displacement = shape->angular_velocity[2] * delta_t;

    rotate_around(shape, x_angular_displacement, y_angular_displacement, z_angular_displacement,
                shape->origin[0], shape->origin[1], shape->origin[2], handle);
}

// find thetas required to rotate vec1 onto vec2
void angle_between_vectors(double *vec1, double *vec2, double *result) {
    result[0] = 0, result[1] = 0, result[2] = 0;

    if (vec1[0] != vec2[0]) { // angle in x/z plane
        double vec1_x_theta = 0;
        if (!(vec1[0] == 0 && vec1[2] == 0)) {
            vec1_x_theta = vec1[0] == 0 ? (vec1[2] < 0 ? -PI / 2.0 : PI / 2.0) : atan(vec1[2] / vec1[0]);
        }
        double vec2_x_theta = 0;
        if (!(vec2[0] == 0 && vec2[2] == 0)) {
            vec2_x_theta = vec2[0] == 0 ? (vec2[2] < 0 ? -PI / 2.0 : PI / 2.0) : atan(vec2[2] / vec2[0]);
        }

        // move to quadrants 2/3 if x negative
        vec1_x_theta += vec1[0] < 0 ? PI : 0;
        vec2_x_theta += vec2[0] < 0 ? PI : 0;

        // use only positive angles
        vec1_x_theta += vec1_x_theta < 0 ? 2 * PI : 0;
        vec2_x_theta += vec2_x_theta < 0 ? 2 * PI : 0;

        // issues with 0 angle b/c -0 is treated as negative
        if (!((vec1_x_theta == 0 || vec1_x_theta == PI) && (vec2_x_theta == 0 || vec2_x_theta == PI))) {
            result[0] = vec2_x_theta - vec1_x_theta;
        }
    }
    if (vec1[1] != vec2[1]) { // angle in x/y plane
        double vec1_y_theta = 0;
        if (!(vec1[0] == 0 && vec1[1] == 0)) {
            vec1_y_theta = vec1[0] == 0 ? (vec1[1] < 0 ? -PI / 2.0 : PI / 2.0) : atan(vec1[1] / vec1[0]);
        }
        double vec2_y_theta = 0;
        if (!(vec2[0] == 0 && vec2[1] == 0)) {
            vec2_y_theta = vec2[0] == 0 ? (vec2[1] < 0 ? -PI / 2.0 : PI / 2.0) : atan(vec2[1] / vec2[0]);
        }

        // move to quadrants 2/3 if x negative
        vec1_y_theta += vec1[0] < 0 ? PI : 0;
        vec2_y_theta += vec2[0] < 0 ? PI : 0;

        // use only positive angles
        vec1_y_theta += vec1_y_theta < 0 ? 2 * PI : 0;
        vec2_y_theta += vec2_y_theta < 0 ? 2 * PI : 0;

        // issues with 0 angle b/c -0 is treated as negative
        if (!((vec1_y_theta == 0 || vec1_y_theta == PI) && (vec2_y_theta == 0 || vec2_y_theta == PI))) {
            result[1] = vec2_y_theta - vec1_y_theta;
        }
    }
    if (vec1[2] != vec2[2]) { // angle in y/z plane
        double vec1_z_theta = 0;
        if (!(vec1[1] == 0 && vec1[2] == 0)) {
            vec1_z_theta = vec1[2] == 0 ? (vec1[1] < 0 ? -PI / 2.0 : PI / 2.0) : atan(vec1[1] / vec1[2]);
        }
        double vec2_z_theta = 0;
        if (!(vec2[1] == 0 && vec2[2] == 0)) {
            vec2_z_theta = vec2[2] == 0 ? (vec2[1] < 0 ? -PI / 2.0 : PI / 2.0) : atan(vec2[1] / vec2[2]);
        }

        // move to quadrants 2/3 if x negative
        vec1_z_theta += vec1[2] < 0 ? PI : 0;
        vec2_z_theta += vec2[2] < 0 ? PI : 0;

        // use only positive angles
        vec1_z_theta += vec1_z_theta < 0 ? 2 * PI : 0;
        vec2_z_theta += vec2_z_theta < 0 ? 2 * PI : 0;

        // issues with 0 angle b/c -0 is treated as negative
        if (!((vec1_z_theta == 0 || vec1_z_theta == PI) && (vec2_z_theta == 0 || vec2_z_theta == PI))) {
            result[2] = vec2_z_theta - vec1_z_theta;
        }
    }
}

void convert_angular_to_linear(Shape *shape, double *shape1_vel_angular, double *linear_vel, double *collision_coords) {
    // find direction; tangent of y-axis to collision point
    double shape1_to_collision[3];
    double s_1_collision_sum = 0;
    for (unsigned int i = 0; i < 3; i++) {
        shape1_to_collision[i] = collision_coords[i] - shape->origin[i];
        s_1_collision_sum += abs(shape1_to_collision[i]);
    }
    double x_sum = abs(shape1_to_collision[0]) + abs(shape1_to_collision[2]) == 0 ? 1 : abs(shape1_to_collision[0]) + abs(shape1_to_collision[2]);
    double s1_angular_x[3] = { 0, 0, 0 };
    s1_angular_x[0] = abs(shape1_to_collision[2]) / x_sum * (shape1_vel_angular[0] >= 0 ? 1 : -1);
    s1_angular_x[2] = abs(shape1_to_collision[0]) / x_sum * (shape1_vel_angular[0] >= 0 ? -1 : 1);
    // angular to linear
    double x_radius = sqrt(pow(s1_angular_x[0] * x_sum, 2) + pow(s1_angular_x[2] * x_sum, 2)); //radius at height
    s1_angular_x[0] *= x_radius * abs(shape1_vel_angular[0]);
    s1_angular_x[2] *= x_radius * abs(shape1_vel_angular[0]);
    linear_vel[0] += s1_angular_x[0];
    linear_vel[2] += s1_angular_x[2];

    // y angular to linear
    double y_sum = abs(shape1_to_collision[0]) + abs(shape1_to_collision[1]) == 0 ? 1 : abs(shape1_to_collision[0]) + abs(shape1_to_collision[1]);
    double s1_angular_y[3] = { 0, 0, 0 };
    s1_angular_y[0] = abs(shape1_to_collision[1]) / y_sum * (shape1_vel_angular[1] >= 0 ? -1 : 1);
    s1_angular_y[1] = abs(shape1_to_collision[0]) / y_sum * (shape1_vel_angular[1] >= 0 ? 1 : -1);
    // angular to linear
    double y_radius = sqrt(pow(s1_angular_y[0] * y_sum, 2) + pow(s1_angular_y[1] * y_sum, 2));
    s1_angular_y[0] *= y_radius * abs(shape1_vel_angular[1]);
    s1_angular_y[1] *= y_radius * abs(shape1_vel_angular[1]);
    linear_vel[0] += s1_angular_y[0];
    linear_vel[1] += s1_angular_y[1];

    // z angular to linear
    double z_sum = abs(shape1_to_collision[1]) + abs(shape1_to_collision[2]) == 0 ? 1 : abs(shape1_to_collision[1]) + abs(shape1_to_collision[2]);
    double s1_angular_z[3] = { 0, 0, 0 };
    s1_angular_z[1] = abs(shape1_to_collision[2]) / z_sum * (shape1_vel_angular[2] >= 0 ? -1 : 1);
    s1_angular_z[2] = abs(shape1_to_collision[1]) / z_sum * (shape1_vel_angular[2] >= 0 ? 1 : -1);
    double z_radius = sqrt(pow(s1_angular_z[1] * z_sum, 2) + pow(s1_angular_z[2] * z_sum, 2));
    s1_angular_z[1] *= z_radius * abs(shape1_vel_angular[2]);
    s1_angular_z[2] *= z_radius * abs(shape1_vel_angular[2]);
    linear_vel[1] += s1_angular_z[1];
    linear_vel[2] += s1_angular_z[2];
}

void convert_linear_to_angular(Shape *shape, double *shape1_linear, double *angular_vel, double *collision_coords) {
    double linear_unit[3] = { shape1_linear[0], shape1_linear[1], shape1_linear[2] };
    double linear_sum = 0;
    for (unsigned int i = 0; i < 3; i++) {
        linear_sum += abs(shape1_linear[i]);
    }
    for (unsigned int i = 0; i < 3; i++) {
        linear_unit[i] /= linear_sum == 0 ? 1 : linear_sum;
    }
    double linear_mag = sqrt(pow(shape1_linear[0], 2) + pow(shape1_linear[1], 2) + pow(shape1_linear[2], 2));

    double shape1_to_collision[3];
    double s_1_collision_sum = 0;
    for (unsigned int i = 0; i < 3; i++) {
        shape1_to_collision[i] = collision_coords[i] - shape->origin[i];
        s_1_collision_sum += abs(shape1_to_collision[i]);
    }

    double x_direction;
    if ((collision_coords[2] >= 0 && shape1_linear[0] >= 0) || (collision_coords[2] < 0 && shape1_linear[0] < 0)) {
        x_direction = 1;
    } else {
        x_direction = -1;
    }
    double x_sum = abs(shape1_to_collision[0]) + abs(shape1_to_collision[2]) == 0 ? 1 : abs(shape1_to_collision[0]) + abs(shape1_to_collision[2]);
    double s1_angular_x[3] = { 0, 0, 0 };
    s1_angular_x[0] = abs(shape1_to_collision[2]) / x_sum * (shape1_linear[0] >= 0 ? 1 : -1);
    s1_angular_x[2] = abs(shape1_to_collision[0]) / x_sum * (shape1_linear[0] >= 0 ? -1 : 1);
    double overlap = dot(linear_unit, s1_angular_x);
    // angular to linear
    double x_radius = sqrt(pow(shape1_to_collision[0], 2) + pow(shape1_to_collision[2], 2)); //radius at height
    angular_vel[0] += x_radius == 0 ? 0 : linear_mag * overlap / (2 * PI * x_radius) * x_direction;

    // y angular to linear
    double y_direction;
    if ((collision_coords[1] >= 0 && shape1_linear[0] < 0) || (collision_coords[1] < 0 && shape1_linear[0] >= 0)) {
        y_direction = 1;
    } else {
        y_direction = -1;
    }
    double y_sum = abs(shape1_to_collision[0]) + abs(shape1_to_collision[1]) == 0 ? 1 : abs(shape1_to_collision[0]) + abs(shape1_to_collision[1]);
    double s1_angular_y[3] = { 0, 0, 0 };
    s1_angular_y[0] = abs(shape1_to_collision[1]) / y_sum * (shape1_linear[1] >= 0 ? -1 : 1);
    s1_angular_y[1] = abs(shape1_to_collision[0]) / y_sum * (shape1_linear[1] >= 0 ? 1 : -1);
    overlap = dot(linear_unit, s1_angular_y);
    // angular to linear
    double y_radius = sqrt(pow(shape1_to_collision[0], 2) + pow(shape1_to_collision[1], 2)); //radius at height
    angular_vel[1] += y_radius == 0 ? 0 : linear_mag * overlap / (2 * PI * y_radius) * y_direction;

    // z angular to linear
    double z_direction;
    if ((collision_coords[2] >= 0 && shape1_linear[1] <= 0) || (collision_coords[2] < 0 && shape1_linear[1] > 0)) {
        z_direction = 1;
    } else {
        z_direction = -1;
    }
    double z_sum = abs(shape1_to_collision[1]) + abs(shape1_to_collision[2]) == 0 ? 1 : abs(shape1_to_collision[1]) + abs(shape1_to_collision[2]);
    double s1_angular_z[3] = { 0, 0, 0 };
    s1_angular_z[1] = abs(shape1_to_collision[2]) / z_sum * (shape1_linear[2] >= 0 ? -1 : 1);
    s1_angular_z[2] = abs(shape1_to_collision[1]) / z_sum * (shape1_linear[2] >= 0 ? 1 : -1);
    overlap = dot(linear_unit, s1_angular_z);
    // angular to linear
    double z_radius = sqrt(pow(shape1_to_collision[1], 2) + pow(shape1_to_collision[2], 2)); //radius at height
    angular_vel[2] += z_radius == 0 ? 0 : linear_mag * overlap / (2 * PI * z_radius) * z_direction;
}

void calculate_hit(Shape *shape1, Shape *shape2, double *collision_coords,
                    double *shape1_vel_input, double *shape2_vel_input, double *shape1_vel_output,
                    double *shape2_vel_output, double *shape1_vel_output_angular, double *shape2_vel_output_angular) {

    // calculate magnitudes
    double origins_vec[] = { shape2->origin[0] - shape1->origin[0], shape2->origin[1] - shape1->origin[1],
                                shape2->origin[2] - shape1->origin[2] };
    double origin_mag = sqrt(pow(origins_vec[0], 2) + pow(origins_vec[1], 2) + pow(origins_vec[2], 2));
    double shape1_collision_mag = sqrt(pow(collision_coords[0] - shape1->origin[0], 2)
                                             + pow(collision_coords[1] - shape1->origin[1], 2)
                                             + pow(collision_coords[2] - shape1->origin[2], 2));
    double shape2_collision_mag = sqrt(pow(collision_coords[0] - shape2->origin[0], 2)
                                        + pow(collision_coords[1] - shape2->origin[1], 2)
                                        + pow(collision_coords[2] - shape2->origin[2], 2));

    // calculate unit vectors
    double ref_x_vec[] = { 1, 0, 0 } ;
    double origins_contact_sum = abs(origins_vec[0]) + abs(origins_vec[1]) + abs(origins_vec[2]);
    double shape1_vel_sum = abs(shape1_vel_input[0]) + abs(shape1_vel_input[1]) + abs(shape1_vel_input[2]);
    double shape2_vel_sum = abs(shape2_vel_input[0]) + abs(shape2_vel_input[1]) + abs(shape2_vel_input[2]);
    double shape1_unit[3];
    double shape2_unit[3];
    double shape1_to_collision[3];
    double shape2_to_collision[3];
    double shape1_ratio = shape1_collision_mag / origin_mag;
    double shape2_ratio = -1 * shape2_collision_mag / origin_mag;
    double s1_collision_sum = 0;
    double s2_collision_sum = 0;
    for (unsigned int i = 0; i < 3; i++) {
        shape1_unit[i] = shape1_vel_sum == 0 ? 0 : shape1_vel_input[i] / shape1_vel_sum;
        shape2_unit[i] = shape2_vel_sum == 0 ? 0 : shape2_vel_input[i] / shape2_vel_sum;

        shape1_to_collision[i] = origins_vec[i] * shape1_ratio;
        s1_collision_sum += abs(shape1_to_collision[i]);
        shape2_to_collision[i] = origins_vec[i] * shape2_ratio;
        s2_collision_sum += abs(shape2_to_collision[i]);

        origins_vec[i] /= origins_contact_sum == 0 ? 1 : origins_contact_sum;
    }
    for (unsigned int i = 0; i < 3; i++) {
        shape1_to_collision[i] /= s1_collision_sum == 0 ? 1 : s1_collision_sum;
        shape2_to_collision[i] /= s2_collision_sum == 0 ? 1 : s2_collision_sum;
    }

    // find angles between unit vectors
    double s1_point_thetas[3];
    angle_between_vectors(shape1_to_collision, ref_x_vec, s1_point_thetas);
    double s2_point_thetas[3];
    angle_between_vectors(shape2_to_collision, ref_x_vec, s2_point_thetas);
    double s1_thetas[3];
    angle_between_vectors(shape1_unit, origins_vec, s1_thetas);
    double s2_thetas[3];
    angle_between_vectors(shape2_unit, origins_vec, s2_thetas);

    // correct tangent vector is the one that points toward the velocity vector
    double s1_thetas_tangent[3] = { 0, 0, 0 };
    for (unsigned int i = 0; i < 3; i++) {
        if (s1_thetas[i] > 0) {
            s1_thetas_tangent[i] = -PI / 2.0;
        } else if (s1_thetas[i] < 0) {
            s1_thetas_tangent[i] = PI / 2.0;
        }
    }
    if (s1_thetas[0] == 0 && s1_thetas[1] == 0 && s1_thetas[2] == 0) {
        s1_thetas_tangent[0] = PI / 2.0;
        s1_thetas_tangent[1] = PI / 2.0;
        s1_thetas_tangent[2] = PI / 2.0;
    }
    double tangent_vec[3];
    rotate_vector(origins_vec, tangent_vec, s1_thetas_tangent[0], s1_thetas_tangent[1], s1_thetas_tangent[2]);

    // calculate dot products for scaling
    double shape_1_aligned_dot = dot(shape1_unit, origins_vec);
    double shape_1_tangent_dot = dot(shape1_unit, tangent_vec);
    double shape_2_aligned_dot = dot(shape2_unit, origins_vec);
    double shape_2_tangent_dot = dot(shape2_unit, tangent_vec);

    // rotate velocity vectors to origin and tangent components
    double shape1_aligned_vector[3] = {shape1_vel_input[0], shape1_vel_input[1], shape1_vel_input[2]};
    rotate_vector(shape1_aligned_vector, shape1_aligned_vector, s1_thetas[0], s1_thetas[1], s1_thetas[2]);
    double shape2_aligned_vector[3] = {shape2_vel_input[0], shape2_vel_input[1], shape2_vel_input[2]};
    rotate_vector(shape2_aligned_vector, shape2_aligned_vector, s1_thetas[0], s1_thetas[1], s1_thetas[2]);
    double shape1_tangent_vector[3] = {shape1_vel_input[0], shape1_vel_input[1], shape1_vel_input[2]};
    rotate_vector(shape1_tangent_vector, shape1_tangent_vector, s1_thetas[0] + s1_thetas_tangent[0],
                s1_thetas[1] + s1_thetas_tangent[1], s1_thetas[2] + s1_thetas_tangent[2]);
    double shape2_tangent_vector[3] = {shape2_vel_input[0], shape2_vel_input[1], shape2_vel_input[2]};
    rotate_vector(shape2_tangent_vector, shape2_tangent_vector, -(s1_thetas[0] + s1_thetas_tangent[0]),
        -(s1_thetas[1] + s1_thetas_tangent[1]), -(s1_thetas[2] + s1_thetas_tangent[2]));

    // scale vectors
    for (unsigned int i = 0; i < 3; i++) {
        shape1_aligned_vector[i] *= abs(shape_1_aligned_dot);
        shape2_aligned_vector[i] *= abs(shape_2_aligned_dot);
        shape1_tangent_vector[i] *= abs(shape_1_tangent_dot);
        shape2_tangent_vector[i] *= abs(shape_2_tangent_dot);
    }

    // angular velocity calculated in terms of linear velocity by using positive x-axis as reference point
    double shape1_normalized_tangent[3] = {shape1_vel_input[0], shape1_vel_input[1], shape1_vel_input[2]};
    double shape2_normalized_tangent[3] = {shape2_vel_input[0], shape2_vel_input[1], shape2_vel_input[2]};
    rotate_vector(shape1_tangent_vector, shape1_normalized_tangent, s1_point_thetas[0], s1_point_thetas[1], s1_point_thetas[2]);
    rotate_vector(shape2_tangent_vector, shape2_normalized_tangent, -s2_point_thetas[0], -s2_point_thetas[1], -s2_point_thetas[2]);

    // calculate velocities using conservation of momentum
    double shape1_lin_angular[3];
    double shape2_lin_angular[3];
    for (unsigned int i = 0; i < 3; i++) {
        shape1_vel_output[i] += (shape1->mass - shape2->mass) / (shape1->mass + shape2->mass) * shape1_aligned_vector[i]
                                    + 2 * shape2->mass / (shape1->mass + shape2->mass) * shape2_aligned_vector[i];
        shape2_vel_output[i] += 2 * shape1->mass / (shape1->mass + shape2->mass) * shape1_aligned_vector[i]
                                    + (shape2->mass - shape1->mass) / (shape1->mass + shape2->mass) * shape2_aligned_vector[i];

        shape1_vel_output[i] += shape1_tangent_vector[i];
        shape2_vel_output[i] += shape2_tangent_vector[i];

        shape1_lin_angular[i] = (shape1->mass - shape2->mass) / (shape1->mass + shape2->mass) * shape1_normalized_tangent[i]
        + 2 * shape2->mass / (shape1->mass + shape2->mass) * shape2_normalized_tangent[i];
        shape2_lin_angular[i] = 2 * shape1->mass / (shape1->mass + shape2->mass) * shape1_normalized_tangent[i]
        + (shape2->mass - shape1->mass) / (shape1->mass + shape2->mass) * shape2_normalized_tangent[i];
    }

    // using positive x-axis as reference, x velocity corresponds to z rotation
    shape1_vel_output_angular[0] += shape1_lin_angular[2] / shape1_collision_mag;
    shape2_vel_output_angular[0] += shape2_lin_angular[2] / shape2_collision_mag * -1;

    shape1_vel_output_angular[1] += shape1_lin_angular[1] / shape1_collision_mag;
    shape2_vel_output_angular[1] += shape2_lin_angular[1] / shape2_collision_mag * -1;

    // using positive x-axis as reference, z velocity corresponds to x rotation
    shape1_vel_output_angular[2] += shape1_lin_angular[0] / shape1_collision_mag;
    shape2_vel_output_angular[2] += shape2_lin_angular[0] / shape2_collision_mag * -1;
}

void update_velocities(Shape *shape1, Shape *shape2, double *collision_coords) {
    // calculate linear velocity contribution
    double shape1_vel[3] = {shape1->velocity[0], shape1->velocity[1], shape1->velocity[2]};
    double shape2_vel[3] = {shape2->velocity[0], shape2->velocity[1], shape2->velocity[2]};
    double shape1_vel_angular[3] = { shape1->angular_velocity[0], shape1->angular_velocity[1], shape1->angular_velocity[2] };
    double shape2_vel_angular[3] = { shape2->angular_velocity[0], shape2->angular_velocity[1], shape2->angular_velocity[2] };

    double shape_1_vel_output[3] = {0, 0, 0};
    double shape_2_vel_output[3] = {0, 0, 0};
    double shape_1_vel_output_angular[3] = {0, 0, 0};
    double shape_2_vel_output_angular[3] = {0, 0, 0};

    calculate_hit(shape1, shape2, collision_coords, shape1_vel, shape2_vel, shape_1_vel_output,
        shape_2_vel_output, shape_1_vel_output_angular, shape_2_vel_output_angular);

    double s1_vel[] = { 0, 0, 0 };
    double s2_vel[] = { 0, 0, 0 };
    convert_angular_to_linear(shape1, shape1_vel_angular, s1_vel, collision_coords);
    convert_angular_to_linear(shape2, shape2_vel_angular, s2_vel, collision_coords);

    double dummy_output[3] = {0.0, 0.0, 0.0}; // placeholder; need to prevent shape's angular velocity from turning into linear
    double dummy_output2[3] = {0.0, 0.0, 0.0};
    calculate_hit(shape1, shape2, collision_coords, s1_vel, s2_vel, dummy_output,
        dummy_output2, shape_1_vel_output_angular, shape_2_vel_output_angular);

    double s1_vel_output[] = { 0, 0, 0};
    double s2_vel_output[] = { 0, 0, 0};
    convert_linear_to_angular(shape1, dummy_output, s1_vel_output, collision_coords);
    convert_linear_to_angular(shape2, dummy_output2, s2_vel_output, collision_coords);
    for (unsigned int i = 0; i < 3; i++) {
        shape_1_vel_output_angular[i] += s1_vel_output[i];
        shape_2_vel_output_angular[i] += s2_vel_output[i];
    }

    // copy results to shapes
    for (unsigned int i = 0; i < 3; i++) {
        shape1->velocity[i] = shape_1_vel_output[i];
        shape2->velocity[i] = shape_2_vel_output[i];
        shape1->angular_velocity[i] = shape_1_vel_output_angular[i];
        shape2->angular_velocity[i] = shape_2_vel_output_angular[i];
    }
}

int triangle_intersection(Shape *shape, unsigned int tri_index, Shape *shape2, unsigned int tri_index2, double *tri1_origin, double *tri2_origin, double *coords)
{
    int intersection_found = 0;
    unsigned int intersection_count = 0;
    unsigned int shape1_side1 = shape->triangles[tri_index].side_vector_index;

    // calculate shape_2 triangle plane
    unsigned int shape2_side1 = shape2->triangles[tri_index2].side_vector_index;
    double normal_vector[3];
    cross(normal_vector,
        shape2->vector_array[shape2_side1][0], shape2->vector_array[shape2_side1][1], shape2->vector_array[shape2_side1][2],
        shape2->vector_array[shape2_side1 + 1][0], shape2->vector_array[shape2_side1 + 1][1], shape2->vector_array[shape2_side1 + 1][2]);

    double current_origin[3] = { tri1_origin[0], tri1_origin[1], tri1_origin[2] };
    for (unsigned int tri_outer = 0; tri_outer < 3; tri_outer++) {
        // find point where triangle side intersects shape2 face
        double t_val = -1 * current_origin[0] * normal_vector[0] + normal_vector[0] * tri2_origin[0] - current_origin[1] * normal_vector[1]
                        + normal_vector[1] * tri2_origin[1] - current_origin[2] * normal_vector[2] + normal_vector[2] * tri2_origin[2];
        double t_val_divisor = normal_vector[0] * shape->vector_array[shape1_side1 + tri_outer][0] +
                                normal_vector[1] * shape->vector_array[shape1_side1 + tri_outer][1] +
                                normal_vector[2] * shape->vector_array[shape1_side1 + tri_outer][2];
        if (t_val_divisor == 0) {
            continue;
        }
        t_val /= t_val_divisor;
        double x_point = t_val * shape->vector_array[shape1_side1 + tri_outer][0] + current_origin[0];
        double y_point = t_val * shape->vector_array[shape1_side1 + tri_outer][1] + current_origin[1];
        double z_point = t_val * shape->vector_array[shape1_side1 + tri_outer][2] + current_origin[2];

        // verify intersection point is in both triangles
        bool point_in_s1 = plane_point_in_triangle(shape->vector_array[shape1_side1], shape->vector_array[shape1_side1 + 1],
                                                    shape->vector_array[shape1_side1 + 2], x_point, y_point, z_point, tri1_origin);
        bool point_in_s2 = plane_point_in_triangle(shape2->vector_array[shape2_side1], shape2->vector_array[shape2_side1 + 1],
                                                    shape2->vector_array[shape2_side1 + 2], x_point, y_point, z_point, tri2_origin);
        if (point_in_s2 && point_in_s1) {
            if (intersection_found <= 1) {
                intersection_found = 1;
            }

            coords[0] = x_point;
            coords[1] = y_point;
            coords[2] = z_point;

            // check if shape1 side is parallel to shape2 face
            double x_mult = normal_vector[0] == 0 ? 0 : abs(shape->vector_array[shape1_side1 + tri_outer][0] / normal_vector[0]);
            double y_mult = normal_vector[1] == 0 ? 0 : abs(shape->vector_array[shape1_side1 + tri_outer][1] / normal_vector[1]);
            double z_mult = normal_vector[2] == 0 ? 0 : abs(shape->vector_array[shape1_side1 + tri_outer][2] / normal_vector[2]);
            bool xy_pass = x_mult == 0 || y_mult == 0 || x_mult == y_mult;
            bool xz_pass = x_mult == 0 || z_mult == 0 || x_mult == z_mult;
            bool yz_pass = y_mult == 0 || z_mult == 0 || y_mult == z_mult;
            if (xy_pass && yz_pass && xz_pass) {
                if (intersection_found <= 2) {
                    intersection_found = 2;
                }
                intersection_count++;

                // parallel vector means collision is line; use midpoint as collision point
                coords[3] = current_origin[0] + shape->vector_array[shape1_side1 + tri_outer][0] / 2.0;
                coords[4] = current_origin[1] + shape->vector_array[shape1_side1 + tri_outer][1] / 2.0;
                coords[5] = current_origin[2] + shape->vector_array[shape1_side1 + tri_outer][2] / 2.0;
            }
        }
        current_origin[0] += shape->vector_array[shape1_side1 + tri_outer][0];
        current_origin[1] += shape->vector_array[shape1_side1 + tri_outer][1];
        current_origin[2] += shape->vector_array[shape1_side1 + tri_outer][2];
    }

    if (intersection_count == 3) { // 3 intersections means plane is parallel
        intersection_found = 3;
        // use average of vertices as collision point
        coords[6] = 0;
        coords[7] = 0;
        coords[8] = 0;
        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
                coords[6 + j] += (current_origin[j] + shape->vector_array[shape1_side1 + i][j]) / 3.0;
                current_origin[j] += shape->vector_array[shape1_side1 + i][j];
            }
        }
    }
    return intersection_found;
}

int calculate_plane_lines(Shape *shape, unsigned int tri_index, Shape *shape2, unsigned int tri_index2, double *coords)
{
    unsigned int shape1_origin_index = shape->triangles[tri_index].origin_vector_index;
    unsigned int shape2_origin_index = shape2->triangles[tri_index2].origin_vector_index;

    double shape1_origin[3];
    double shape2_origin[3];
    for (unsigned int i = 0; i < 3; i++) {
        shape1_origin[i] = shape->origin[i] + shape->vector_array[shape1_origin_index][i];
        shape2_origin[i] = shape2->origin[i] + shape2->vector_array[shape2_origin_index][i];
    }

    int shape1_pass = triangle_intersection(shape, tri_index, shape2, tri_index2, shape1_origin, shape2_origin, coords);
    int shape2_pass = triangle_intersection(shape2, tri_index2, shape, tri_index, shape2_origin, shape1_origin, coords);

    if (shape1_pass > 0 || shape2_pass > 0) {
        return max(shape1_pass, shape2_pass);
    }
    return 0;
}

// functor for thrust
struct tri_collision_check {
    Shape *shape1;
    Shape *shape2;
    double *collision_coords;

    tri_collision_check(Shape *shape1, Shape *shape2, double *collision_coords) : shape1(shape1), shape2(shape2), collision_coords(collision_coords) {}

    int operator()(const int& input1, const int& input2) {
        unsigned int tri_index1 = input1 / 12;
        unsigned int tri_index2 = input2 % 12;
        int collision = calculate_plane_lines(shape1, tri_index1, shape2, tri_index2, collision_coords);
        return collision;
    }
};

bool check_collision(Shape *shape1, Shape *shape2)
{
    // find shape bounding box
    double shape_1_x_min_max[2] = {shape1->origin[0], shape1->origin[0]};
    double shape_1_y_min_max[2] = {shape1->origin[1], shape1->origin[1]};
    double shape_1_z_min_max[2] = {shape1->origin[2], shape1->origin[2]};
    shape_min_max_coords(shape1, shape_1_x_min_max, shape_1_y_min_max, shape_1_z_min_max);

    double shape_2_x_min_max[2] = {shape2->origin[0], shape2->origin[0]};
    double shape_2_y_min_max[2] = {shape2->origin[1], shape2->origin[1]};
    double shape_2_z_min_max[2] = {shape2->origin[2], shape2->origin[2]};
    shape_min_max_coords(shape2, shape_2_x_min_max, shape_2_y_min_max, shape_2_z_min_max);

    bool impossible_x_intersection = (shape_1_x_min_max[1] < shape_2_x_min_max[0] || shape_1_x_min_max[0] > shape_2_x_min_max[1]);
    bool impossible_y_intersection = (shape_1_y_min_max[1] < shape_2_y_min_max[0] || shape_1_y_min_max[0] > shape_2_y_min_max[1]);
    bool impossible_z_intersection = (shape_1_z_min_max[1] < shape_2_z_min_max[0] || shape_1_z_min_max[0] > shape_2_z_min_max[1]);
    if (impossible_x_intersection || impossible_y_intersection || impossible_z_intersection) {
        return false;
    }

    thrust::host_vector<int> shape_1_vec(shape1->triangle_count * shape2->triangle_count);
    thrust::host_vector<int> shape_2_vec(shape1->triangle_count * shape2->triangle_count);
    thrust::host_vector<int> output_vec(shape1->triangle_count * shape2->triangle_count);

    thrust::sequence(thrust::host, shape_1_vec.begin(), shape_1_vec.end());
    thrust::sequence(thrust::host, shape_2_vec.begin(), shape_2_vec.end());

    double collision_coords[9];
    thrust::transform(thrust::host, shape_1_vec.begin(), shape_1_vec.end(), shape_2_vec.begin(),
                        output_vec.begin(), tri_collision_check(shape1, shape2, collision_coords));

    int collision = *thrust::max_element(output_vec.begin(), output_vec.end());
    if (collision > 0) {
        double passed_collision[3];
        for (unsigned int i = 0; i < 3; i++) {
            passed_collision[i] = collision_coords[(collision - 1) * 3 + i];
        }

        update_velocities(shape1, shape2, passed_collision);
        printf("collision detected\n");

        // shift_vector(shape1->origin, -200, 0, 0);

        return true;
    }
    return false;
}

void model_movement(Shape **input_shapes, unsigned int shape_count, int x_rotate, int y_rotate, int z_rotate,
                    double initial_timestep, double final_timestep) {

    // rotate around origin
    double axis_x = input_shapes[0]->origin[0];
    double axis_y = input_shapes[0]->origin[1];
    double axis_z = input_shapes[0]->origin[2];

    cublasHandle_t cnpHandle;
    cublasStatus_t status = cublasCreate_v2(&cnpHandle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("Cublas error: %d\n", status);
    }
    rotate_around(input_shapes[0], degrees_to_radians(x_rotate), degrees_to_radians(y_rotate), degrees_to_radians(z_rotate),
                    axis_x, axis_y, axis_z, &cnpHandle);

    for (unsigned int i = 0; i < shape_count; i++) {
        adjust_movement(input_shapes[i], initial_timestep, final_timestep);
        adjust_rotation(input_shapes[i], initial_timestep, final_timestep, &cnpHandle);
    }

    check_collision(input_shapes[0], input_shapes[1]);

    cublasDestroy_v2(cnpHandle);
}

void draw_shapes(Pixel *bitmap, Shape **input_shapes, unsigned int shape_count) {
    double side_len = 10; // placeholder

    int texture_size = (int) (side_len * side_len);
    Pixel *blue_bitmap = new Pixel[texture_size];
    for (int i = 0; i < texture_size; i++) {
        blue_bitmap[i].blue = 255;
        blue_bitmap[i].green = 0;
        blue_bitmap[i].red = 0;
    }
    Pixel *green_bitmap = new Pixel[texture_size];
    for (int i = 0; i < texture_size; i++) {
        green_bitmap[i].blue = 0;
        green_bitmap[i].green = 255;
        green_bitmap[i].red = 0;
    }
    Pixel *red_bitmap = new Pixel[texture_size];
    for (int i = 0; i < texture_size; i++) {
        red_bitmap[i].blue = 0;
        red_bitmap[i].green = 0;
        red_bitmap[i].red = 255;
    }
    Pixel *pink_bitmap = new Pixel[texture_size];
    for (int i = 0; i < texture_size; i++) {
        pink_bitmap[i].blue = 255;
        pink_bitmap[i].green = 0;
        pink_bitmap[i].red = 255;
    }
    Pixel *yellow_bitmap = new Pixel[texture_size];
    for (int i = 0; i < texture_size; i++) {
        yellow_bitmap[i].blue = 0;
        yellow_bitmap[i].green = 255;
        yellow_bitmap[i].red = 255;
    }
    Pixel *cyan_bitmap = new Pixel[texture_size];
    for (int i = 0; i < texture_size; i++) {
        cyan_bitmap[i].blue = 255;
        cyan_bitmap[i].green = 255;
        cyan_bitmap[i].red = 0;
    }

    double **depth_map = new double*[IMAGE_HEIGHT];
    for (int i = 0; i < IMAGE_HEIGHT; i++) {
        depth_map[i] = new double[IMAGE_WIDTH];
    }
    for (int i = 0; i < (IMAGE_HEIGHT); i++) {
        for (int j = 0; j < (IMAGE_WIDTH); j++) {
            depth_map[i][j] = -1.0;
        }
    }

    // zero passed bitmap
    // allocate GPU memory
    Npp8u *srcPtr;
    Npp8u *dstPtr;
    cudaMalloc(&srcPtr, IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(uint8_t) * 3);
    cudaMalloc(&dstPtr, IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(uint8_t) * 3);

    // copy bitmap to GPU memory
    cudaMemcpy(dstPtr, bitmap, IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(uint8_t) * 3, cudaMemcpyHostToDevice);

    // Use NPP to subtract the base bitmap from itself
    int bytes_per_row = sizeof(uint8_t) * 3 * IMAGE_WIDTH;
    NppiSize roiSize;
    roiSize.width = IMAGE_WIDTH;
    roiSize.height = IMAGE_HEIGHT;
    cudaDeviceSynchronize();

    NppStatus status = nppiSub_8u_C3IRSfs(dstPtr, bytes_per_row, dstPtr, bytes_per_row, roiSize, 0);
    if (status != NPP_SUCCESS) {
        printf("Cublas error: %d\n", status);
    }

    // copy result from GPU memory
    cudaMemcpy(bitmap, dstPtr, IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(uint8_t) * 3, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    // free GPU memory
    cudaFree(dstPtr);
    cudaFree(srcPtr);

    // iterate over input shape array
    for (unsigned int i = 0; i < shape_count; i++) {
        Shape *cube = input_shapes[i];

        // placeholder pending mapping changes for bitmap
        cube->triangles[0].bitmap = green_bitmap;
        cube->triangles[1].bitmap = green_bitmap;
        cube->triangles[2].bitmap = blue_bitmap;
        cube->triangles[3].bitmap = blue_bitmap;
        cube->triangles[4].bitmap = red_bitmap;
        cube->triangles[5].bitmap = red_bitmap;
        cube->triangles[6].bitmap = pink_bitmap;
        cube->triangles[7].bitmap = pink_bitmap;
        cube->triangles[8].bitmap = yellow_bitmap;
        cube->triangles[9].bitmap = yellow_bitmap;
        cube->triangles[10].bitmap = cyan_bitmap;
        cube->triangles[11].bitmap = cyan_bitmap;

        for (unsigned int triangle_index = 0; triangle_index < cube->triangle_count; triangle_index++) {
            draw_shape(cube, triangle_index, bitmap, depth_map);
        }
        for (unsigned int texture_index = 0; texture_index < cube->triangle_count; texture_index++) {
            draw_texture(cube, texture_index, bitmap, depth_map);
        }
    }

    // clean up memory
    delete[] blue_bitmap;
    delete[] green_bitmap;
    delete[] red_bitmap;
    delete[] pink_bitmap;
    delete[] yellow_bitmap;
    delete[] cyan_bitmap;
    for (int i = 0; i < (IMAGE_HEIGHT); i++) {
        delete[] depth_map[i];
    }
    delete[] depth_map;
}

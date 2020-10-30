#include "ray_differentials.h"

// updates the ray differentials so that they match the t value for the intersection position.
// See equations in sections 3.1.1 (transfer) of the quoted paper
void transfer_ray_differentials(Ray &ray, glm::vec3 normal)
{
    glm::vec3 N = glm::normalize(normal);
    glm::vec3 D = glm::normalize(ray.direction);

    float dt_dx = -glm::dot(ray.dP_dx + ray.t * ray.dD_dx, N) / glm::dot(D, N);
    float dt_dy = -glm::dot(ray.dP_dy + ray.t * ray.dD_dy, N) / glm::dot(D, N);

    ray.dP_dx = (ray.dP_dx + ray.t * ray.dD_dx) + dt_dx * D;
    ray.dP_dy = (ray.dP_dy + ray.t * ray.dD_dy) + dt_dy * D;
}

// updates the ray differentials when the ray gets refelected
// See the equations in sections 3.1.2 (reflection) of the quoted paper
// dN_dx, and dN_dy are the derivatives of the point normal in terms of the x and y positions of the screen
void reflect_ray_differentials(Ray &ray, glm::vec3 normal, glm::vec3 direction_before_reflection, glm::vec3 dN_dx, glm::vec3 dN_dy)
{
    glm::vec3 N = glm::normalize(normal);
    glm::vec3 D = glm::normalize(direction_before_reflection);

    float dDN_dx = glm::dot(ray.dD_dx, N) + glm::dot(D, dN_dx);
    float dDN_dy = glm::dot(ray.dD_dy, N) + glm::dot(D, dN_dy);

    ray.dD_dx = ray.dD_dx - 2.0f * (glm::dot(D, N) * dN_dx + dDN_dx * N);
    ray.dD_dy = ray.dD_dy - 2.0f * (glm::dot(D, N) * dN_dy + dDN_dy * N);
}

// c is the reference point for the barycentric coordinate
// a and b form the triange
// p is the point for which we are computing the coordinate
// p is the partial derivative of x in terms of the axis we are interested in
// area is the signed area of the triangle based parallelogram, computed such that the coordinates come out postive
float computeDerivativeOfBarycentricCoordinate(glm::vec3 a, glm::vec3 b, glm::vec3 p, glm::vec3 p_derivative, float area)
{
    glm::vec3 term1 = glm::cross(p_derivative, p - b) + glm::cross(p - a, p_derivative);
    glm::vec3 term2 = glm::cross(a - p, b - p);

    float nominator = glm::dot(term1, term2) + glm::dot(term2, term1);
    float denominator = 2 * area * glm::sqrt(glm::dot(term2, term2));

    return nominator / denominator;
}

// computes the the partial derivatitive of the normal on the intersection point P of the triangle
// if p_derivative is the partial derivative in terms of x, then the derivative of the normal will also be in terms of x
// else the deribative of the normal will be in terms of y
glm::vec3 computeNormalPartialDerivativeInInterpolatedTrianglePoint(Vertex v0, Vertex v1, Vertex v2, glm::vec3 p, glm::vec3 p_derivative)
{
    float area = glm::length(glm::cross(v2.p - v0.p, v1.p - v0.p));

    float alpha_derivative = computeDerivativeOfBarycentricCoordinate(v2.p, v1.p, p, p_derivative, area);
    float beta_derivative = computeDerivativeOfBarycentricCoordinate(v0.p, v2.p, p, p_derivative, area);
    float gamma_derivative = computeDerivativeOfBarycentricCoordinate(v1.p, v0.p, p, p_derivative, area);

    glm::vec3 N_alpha = glm::normalize(v0.n);
    glm::vec3 N_beta = glm::normalize(v1.n);
    glm::vec3 N_gamma = glm::normalize(v2.n);

    return (alpha_derivative * N_alpha) + (beta_derivative * N_beta) + (gamma_derivative * N_gamma);
}

// computes the partial derivative of the texture coordinate on the point p of the triangle
// if p_derivative is the partial derivative in terms of x, then the derivative of the normal will also be in terms of x
// else the deribative of the normal will be in terms of y
glm::vec2 computeTexturePartialDerivativeInInterpolatedTrianglePoint(Vertex v0, Vertex v1, Vertex v2, glm::vec3 p, glm::vec3 p_derivative)
{
    float area = glm::length(glm::cross(v2.p - v0.p, v1.p - v0.p));

    float alpha_derivative = computeDerivativeOfBarycentricCoordinate(v2.p, v1.p, p, p_derivative, area);
    float beta_derivative = computeDerivativeOfBarycentricCoordinate(v0.p, v2.p, p, p_derivative, area);
    float gamma_derivative = computeDerivativeOfBarycentricCoordinate(v1.p, v0.p, p, p_derivative, area);

    glm::vec2 T_alpha = v0.texCoord;
    glm::vec2 T_beta = v1.texCoord;
    glm::vec2 T_gamma = v2.texCoord;

    return (alpha_derivative * T_alpha) + (beta_derivative * T_beta) + (gamma_derivative * T_gamma);
}

// computes the partial derivative of the normal on the point P of the sphere
glm::vec3 computeNormalPartialDerivativeInSphere(float radius, glm::vec3 p_derivative)
{
    return p_derivative / radius;
}

// transfers and reflects the ray differentials when the ray hits some surface
// to be called before actually reflecting the ray
void tranfer_and_reflect_ray_differentials(Ray &ray, HitInfo &hitInfo)
{
    // the object was hit, so transfer the ray differentials to the intersection point
    transfer_ray_differentials(ray, hitInfo.normal);

    // find the partial derivatives of the normals
    glm::vec3 dN_dx, dN_dy;

    if (hitInfo.is_triangle)
    {
        // for interpolated triangles
        dN_dx = computeNormalPartialDerivativeInInterpolatedTrianglePoint(hitInfo.intersected_triangle[0], hitInfo.intersected_triangle[1], hitInfo.intersected_triangle[2], hitInfo.hitPoint, ray.dP_dx);
        dN_dy = computeNormalPartialDerivativeInInterpolatedTrianglePoint(hitInfo.intersected_triangle[0], hitInfo.intersected_triangle[1], hitInfo.intersected_triangle[2], hitInfo.hitPoint, ray.dP_dy);
    }
    else
    {
        // for spheres
        dN_dx = computeNormalPartialDerivativeInSphere(hitInfo.intersected_sphere.radius, ray.dP_dx);
        dN_dy = computeNormalPartialDerivativeInSphere(hitInfo.intersected_sphere.radius, ray.dP_dy);
    }

    reflect_ray_differentials(ray, hitInfo.normal, ray.direction, dN_dx, dN_dy);
}

// compute the ratio between the number of screen pixels and the number of texture pixels for the provided mipmap level
// returns a negative value in the case of an invalid level or uninitialized mipmap, or if the operation is unsuccesful for some other reasons.
float computeLevelOfDetails(const Ray &ray, const HitInfo &hitInfo)
{
    // if the the shape is not a triangle we do not support mipmapping, so return the default value and throw a warning
    if (!hitInfo.is_triangle)
    {
        std::cerr << "Mipmapping not supported on spheres!" << std::endl;
        return 0.0f;
    }

    float deltaX = 1.0f; // corresponds to moving one pixel in the x direction
    float deltaY = 1.0f; // corresponds to moving one pixel in the y direction

    glm::vec2 dT_dx = computeTexturePartialDerivativeInInterpolatedTrianglePoint(hitInfo.intersected_triangle[0], hitInfo.intersected_triangle[1], hitInfo.intersected_triangle[2], hitInfo.hitPoint, ray.dP_dx);
    glm::vec2 dT_dy = computeTexturePartialDerivativeInInterpolatedTrianglePoint(hitInfo.intersected_triangle[0], hitInfo.intersected_triangle[1], hitInfo.intersected_triangle[2], hitInfo.hitPoint, ray.dP_dy);

    glm::vec2 deltaTx = deltaX * dT_dx;
    glm::vec2 deltaTy = deltaY * dT_dy;

    // oversampling is considered  so we are choosing the larger differential
    // in case of negative results, return 0.
    return glm::max(0.0f, glm::log2(glm::max(glm::length(deltaTx), glm::length(deltaTy))));
}

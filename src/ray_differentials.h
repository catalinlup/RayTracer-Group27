

#pragma once
#include "ray.h"
#include "ray_tracing.h"
#include "scene.h"
#include "disable_all_warnings.h"
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/common.hpp>
#include <glm/vec4.hpp>
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <iostream>



// Methods used for Ray Differentials. See paper https://graphics.stanford.edu/papers/trd/trd_jpg.pdf
// The Ray differentials are necessary to select the right level of mipmaps

// transfers and reflects the ray differentials when the ray hits some surface
void tranfer_and_reflect_ray_differentials(Ray &ray, HitInfo &hitInfo);

// updates the ray differentials so that they match the t value for the intersection position.
// See equations in sections 3.1.1 (transfer) of the quoted paper
void transfer_ray_differentials(Ray &ray, glm::vec3 normal);

// updates the ray differentials when the ray gets refelected
// See the equations in sections 3.1.2 (reflection) of the quoted paper
// dN_dx, and dN_dy are the derivatives of the point normal in terms of the x and y positions of the screen
void reflect_ray_differentials(Ray &ray, glm::vec3 normal, glm::vec3 direction_before_reflection, glm::vec3 dN_dx, glm::vec3 dN_dy);

// computes the derivative of the specified barycentric coordinate for p
float computeDerivativeOfBarycentricCoordinate(glm::vec3 a, glm::vec3 b, glm::vec3 p, glm::vec3 p_derivative, float area);

// computes the the partial derivatitive of the normal on the point P of the triangle
// if p_derivative is the partial derivative in terms of x, then the derivative of the normal will also be in terms of x
// else the deribative of the normal will be in terms of y
glm::vec3 computeNormalPartialDerivativeInInterpolatedTrianglePoint(Vertex v0, Vertex v1, Vertex v2, glm::vec3 p, glm::vec3 p_derivative);

// computes the partial derivative of the texture coordinate on the point p of the triangle
// if p_derivative is the partial derivative in terms of x, then the derivative of the normal will also be in terms of x
// else the deribative of the normal will be in terms of y
glm::vec2 computeTexturePartialDerivativeInInterpolatedTrianglePoint(Vertex v0, Vertex v1, Vertex v2, glm::vec3 p, glm::vec3 p_derivative);

// computes the partial derivative of the normal on the point P of the sphere
glm::vec3 computeNormalPartialDerivativeInSphere(float radius, glm::vec3 p_derivative);

// makes use of the ray differentials to compute the required level of details for a mipmap
// the intuition behind the formula is that it computes how many texels we move in the texture space, if we move one pixel in the screen space
// and takes the logarithm of that to get a floating point level
float computeLevelOfDetails(const Ray &ray, const HitInfo &hitInfo);


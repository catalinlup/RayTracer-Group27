#pragma once
#include "disable_all_warnings.h"
// Suppress warnings in third-party code.
DISABLE_WARNINGS_PUSH()
#include <glm/vec3.hpp>
#include <glm/common.hpp>
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <limits>

struct Ray {
    glm::vec3 origin { 0.0f };
    glm::vec3 direction { 0.0f, 0.0f, -1.0f };
    float t { std::numeric_limits<float>::max() };

    // ray differentials used in mipmap selection

    // see equations (8) from this paper https://graphics.stanford.edu/papers/trd/trd_jpg.pdf
    glm::vec3 dD_dx = (glm::dot(direction, direction) * right - glm::dot(direction, right) * direction) / glm::pow(glm::dot(direction, direction), 1.5f);
    glm::vec3 dD_dy = (glm::dot(direction, direction) * up - glm::dot(direction, up) * direction) / glm::pow(glm::dot(direction, direction), 1.5f);
    glm::vec3 dP_dx = glm::vec3(0);
    glm::vec3 dP_dy = glm::vec3(0);

    // defines the right direction of the screen pane
    glm::vec3 right = glm::vec3(1, 0, 0);

    // defines the up direction of the scree pane
    glm::vec3 up = glm::vec3(0, -1, 0);
};

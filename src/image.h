#pragma once
#include "disable_all_warnings.h"
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <vector>

enum class OutOfBoundsRule {
    Border,
    Clamp,
    Repeat
};

enum class TextureFiltering {
    Bilinear,
    NearestNeighbor
};

class Image {
public:
    Image(const std::filesystem::path& filePath);

    glm::vec3 getPixel(const glm::vec2& textureCoordinates) const;

    // sets the border color of the image
    void setBorderColor(const glm::vec3 color);

    // sets the out of bounds rule of the x axis
    void setOutOfBoundsRuleX(const OutOfBoundsRule rule);

    // sets the out of bounds rule of the y axis
    void setOutOfBoundsRuleY(const OutOfBoundsRule rule);

    // used to set the texture filtering method to be used when drawing the texture.
    void setTextureFilteringMethod(const TextureFiltering method);

private:

    // the texture coordinates are in the space [0, 1] x [0, 1], with (0, 0) being the lower-left corner, and (1, 1) the upper right corner
    // this functions maps them to the space [0, m_width] x [0, m_height]
    glm::vec2 toImageCoordinates(glm::vec2 textureCoordinates) const;

    // applies the clamp and repeat procedures on the texture coordinates, on the basis of the selected modes of operations
    glm::vec2 dealWithClampingAndRepeating(glm::vec2 textureCoordinates) const;

    // applies clamping procedure on the provided texture coordinate
    float clampTextureCoordinate(const float coordinate) const;

    // applies the repeating procedure on the provided texture coordinate
    float repeatTextureCoordinate(const float coordinate) const;

    // applies either the clamp or the repeat procedure based on the provided rule
    // if none of the rules were selected, returns the original coordinate
    float clampRepeatTextureCoordinate(const float coordinate, OutOfBoundsRule rule) const;

    // returns true if the texture coordinate is out of bounds,
    // else it returns false
    bool isOutOfBounds(const float coordinate) const;

    // turns the provided (x, y) coordinates into an index that can be used in the m_pixels vector.
    unsigned int linearize2DCoordinates(unsigned int x, unsigned int y) const;

    // turns the provided image coordinates into an index for the pixel array using the nearest neighbor method
    glm::vec3 nearestNeighbor(glm::vec2 imageCoordinates) const;

    // turns the provided image coordinates into an index for the pixel array using the bilinear interpolation method
    glm::vec3 bilinearInterpolation(glm::vec2 imageCoordinates) const;

    // linearly intepolates between 2 colors
    // low - the lower position
    // high - the higher position
    // p - the position for which we are interpolating, p is in [low, high]
    // color_low - the color corresponding to the low position
    // color_high - the color corresponding to the high position
    glm::vec3 linearInterpolation(float low, float high, glm::vec3 color_low, glm::vec3 color_high, float p) const;

    int m_width, m_height;
    std::vector<glm::vec3> m_pixels;

    // the out of bound rule to be applied for the x axis
    // the default is border, i.e. it returns the border color
    OutOfBoundsRule _outOfBoundRuleX = OutOfBoundsRule::Border;

    // the out of bound rule to be applied for the y axis
    // the default is border, i.e. it returns the border color
    OutOfBoundsRule _outOfBoundRuleY = OutOfBoundsRule::Border;

    // the border color to be used when the texture coordinate is out of bounds and the rule is set to Border
    // the default is 0
    glm::vec3 _borderColor = glm::vec3(0);

    // the texture filtering method to be used when rendering textures.
    // the default method is nearest neighbor
    TextureFiltering _filteringMethod = TextureFiltering::NearestNeighbor;
};

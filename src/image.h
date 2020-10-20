#pragma once
#include "disable_all_warnings.h"
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/common.hpp>
#include <glm/vec4.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <vector>
#include <algorithm>

enum class OutOfBoundsRule {
    Border,
    Clamp,
    Repeat
};

enum class TextureFiltering
{
    NearestNeighbor,
    Bilinear,
    MipMappingNearestLevel,
    MipMappingTrilinearFiltering
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

    // FOR MIPMAPPING

    // Returns true if the provided texture supports mipmapping, i.e. if the texture is 2^N x 2^N
    bool canUseMipmapping() const;

    // Returns true if the mipmap is initialized, otherwise false
    bool isMipmapInit() const;

    // returns the number of mipmap levels, -1 if the mipmap is not initialized
    int getNumMipmapLevels() const;

    // returns the currently used level of the mipmap.
    int getCurrMipmapLevel() const;



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
    unsigned int linearize2DCoordinates(unsigned int x, unsigned int y, unsigned int width) const;

    // turns the provided image coordinates into an index for the pixel array using the nearest neighbor method
    glm::vec3 nearestNeighbor(glm::vec2 imageCoordinates, int level, unsigned int width, unsigned int height) const;

    // turns the provided image coordinates into an index for the pixel array using the bilinear interpolation method
    glm::vec3 bilinearInterpolation(glm::vec2 imageCoordinates, int level, unsigned int width, unsigned int height) const;

    // linearly intepolates between 2 colors
    // low - the lower position
    // high - the higher position
    // p - the position for which we are interpolating, p is in [low, high]
    // color_low - the color corresponding to the low position
    // color_high - the color corresponding to the high position
    glm::vec3 linearInterpolation(float low, float high, glm::vec3 color_low, glm::vec3 color_high, float p) const;

    int m_width, m_height;
    std::vector<glm::vec3> m_pixels;

    int curr_level = 0; // the currently selected mimpap level. If no mipmapping is used, it is always set to 0.
    int curr_width, curr_height; // the size of the currently drawn texture. If no mipmapping is used, it is always set m_width, m_height
    std::vector<glm::vec3> curr_texture; // texture that is currently selected for drawing. If no mipmapping is used, it is always set to m_pixels

    bool _mipmap_init = false; // true if the mipmap has been initialized
    std::vector<std::vector<glm::vec3>> _mipmap; // contains all versions of the texture, from the original to the texture of size 1 x 1. If the mipmap was not initialized, it only contains the original texture

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


    // METHODS USED IN MIPMAPPING

    // creates a new texture of reduced resolution corresponding to the next level in the mipmap hierarchy
    // the texture should be 2^n x 2^n
    void getReducedResolutionTexture(int or_width, int or_height, std::vector<glm::vec3>& original, int& re_width, int& re_height, std::vector<glm::vec3>& reduced) const;

    // if the mipmap can be build, and wasn't already initialized, initialize it
    void initMipmap(); 

    // used only for debugging purposes
    void printMipmap() const;

    // switches to another mipmap level. Returns false and does not do anything if it was unsuccesful, i.e. if the mipmap was not initialized or the level is invalid
    bool switchToLevel(int level);


    void printMipmap();
};

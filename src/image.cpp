#include "image.h"
#include "disable_all_warnings.h"
DISABLE_WARNINGS_PUSH()
#include <stb_image.h>
DISABLE_WARNINGS_POP()
#include <cassert>
#include <exception>
#include <iostream>
#include <string>
#include <cmath>

// sets the border color of the image
void Image::setBorderColor(const glm::vec3 color) {
    _borderColor = color;
}

// sets the out of bounds rule of the x axis
void Image::setOutOfBoundsRuleX(const OutOfBoundsRule rule) {
    _outOfBoundRuleX = rule;
}

// sets the out of bounds rule of the y axis
void Image::setOutOfBoundsRuleY(const OutOfBoundsRule rule) {
    _outOfBoundRuleY = rule;
}

// used to set the texture filtering method to be used when drawing the texture.
void Image::setTextureFilteringMethod(const TextureFiltering method) {
    _filteringMethod = method;
}

Image::Image(const std::filesystem::path& filePath)
{
    if (!std::filesystem::exists(filePath)) {
        std::cerr << "Texture file " << filePath << " does not exists!" << std::endl;
        throw std::exception();
    }

    const auto filePathStr = filePath.string(); // Create l-value so c_str() is safe.
    int numChannels;
    stbi_uc* pixels = stbi_load(filePathStr.c_str(), &m_width, &m_height, &numChannels, STBI_rgb);

    if (numChannels < 3) {
        std::cerr << "Only textures with 3 or more color channels are supported. " << filePath << " has " << numChannels << " channels" << std::endl;
        throw std::exception();
    }
    if (!pixels) {
        std::cerr << "Failed to read texture " << filePath << " using stb_image.h" << std::endl;
        throw std::exception();
    }

    std::cout << "Num channels: " << numChannels << std::endl;
    for (size_t i = 0; i < m_width * m_height * numChannels; i += numChannels) {
        m_pixels.emplace_back(pixels[i + 0] / 255.0f, pixels[i + 1] / 255.0f, pixels[i + 2] / 255.0f);
    }

    stbi_image_free(pixels);
}

glm::vec3 Image::getPixel(const glm::vec2& textureCoordinates) const
{
    // deal with the out of bounds cases

    // in case of the border mode:
    if((_outOfBoundRuleX == OutOfBoundsRule::Border) && isOutOfBounds(textureCoordinates.x))
        return _borderColor;
    
    if((_outOfBoundRuleY == OutOfBoundsRule::Border) && isOutOfBounds(textureCoordinates.y))
        return _borderColor;

    // in case of clamping and/or repeating

    glm::vec2 inBoundsTextureCoordinates = dealWithClampingAndRepeating(textureCoordinates);

    glm::vec2 imageCoordinates = toImageCoordinates(inBoundsTextureCoordinates);
    

    // return the color based on one of the 2 filtering methods
    if(_filteringMethod == TextureFiltering::NearestNeighbor)
        return nearestNeighbor(imageCoordinates);
    else if(_filteringMethod == TextureFiltering::Bilinear)
        return bilinearInterpolation(imageCoordinates);

    return glm::vec3(1);
    
}

// the texture coordinates are in the space [0, 1] x [0, 1], with (0, 0) being the lower-left corner, and (1, 1) the upper right corner
// this functions maps them to the space [0, m_width - 1] x [0, m_height - 1]
glm::vec2 Image::toImageCoordinates(glm::vec2 textureCoordinates) const {
    // the y value must be inverted, since the origin is at the top of the screen for an image
    return glm::vec2(textureCoordinates.x * (m_width - 1), (1.0f - textureCoordinates.y) * (m_height - 1));
}



// applies clamping procedure on the provided texture coordinate
float Image::clampTextureCoordinate(const float coordinate) const{
    if (coordinate > 1)
        return 1;
    
    if (coordinate < 0)
        return 0;

    return coordinate;
}

// aplies the repeating procedure on the provided texture coordinate
float Image::repeatTextureCoordinate(const float coordinate) const{
    if(isOutOfBounds(coordinate)) {
        float toRet = std::abs(coordinate);
        float val = toRet - std::floor(toRet);
        // for integers like (2.0f, 3.0f, 4.0f) return 1
        if(val < 1e-7)
            return 1.0f;
        
        return val;
    }

    return coordinate;
}

// returns true if the texture coordinate is out of bounds,
// else it returns false
bool Image::isOutOfBounds(const float coordinate) const {
    if (coordinate < 0 || coordinate > 1)
        return true;

    return false;
}

// applies the clamp and repeat procedures on the texture coordinates, on the basis of the selected modes of operations
glm::vec2 Image::dealWithClampingAndRepeating(glm::vec2 textureCoordinates) const {

    // the x axis
   textureCoordinates.x = clampRepeatTextureCoordinate(textureCoordinates.x, _outOfBoundRuleX);

   // the y axis
   textureCoordinates.y = clampRepeatTextureCoordinate(textureCoordinates.y, _outOfBoundRuleY);

   return textureCoordinates;

}

// applies either the clamp or the repeat procedure based on the provided rule
// if none of the rules were selected, returns the original coordinate
float Image::clampRepeatTextureCoordinate(const float coordinate, OutOfBoundsRule rule) const {
    if (rule == OutOfBoundsRule::Clamp)
        return clampTextureCoordinate(coordinate);

    if (rule == OutOfBoundsRule::Repeat)
        return repeatTextureCoordinate(coordinate);
    
    return coordinate;
}

// turns the provided (x, y) coordinates into an index that can be used in the m_pixels vector.
unsigned int Image::linearize2DCoordinates(unsigned int x, unsigned int y) const {
    // x is the column and y is the row!
    return y * m_width + x;
}   

// returns the corresponding texture color after performing the nearest neighbor approximation
glm::vec3 Image::nearestNeighbor(glm::vec2 imageCoordinates) const {
    unsigned int x = round(imageCoordinates.x);
    unsigned int y = round(imageCoordinates.y);

    // correction, in case some error has occured, to avoid index out of bounds
    if(x >= m_width) {
        x = m_width - 1;
        std::cerr << "Nearest Neighbor: The obtained x was larger than the image width" << std::endl;
    }

    if (x < 0) {
        x = 0;
        std::cerr << "Nearest Neighbor: The obtained x was smaller than 0" << std::endl;
    }

    if (y >= m_height) {
        y = m_height - 1;
        std::cerr << "Nearest Neighbor: The obtained y was larger than the image height" << std::endl;
    }

    if (y < 0) {
        y = 0;
        std::cerr << "Nearest Neighbor: The obtained y was smaller than 0" << std::endl;
    }

    return m_pixels[linearize2DCoordinates(x, y)];
}

// returns the corresponding texture color after applying bilinear interpolation
glm::vec3 Image::bilinearInterpolation(glm::vec2 imageCoordinates) const
{
    float x_low = std::floor(imageCoordinates.x);
    float x_high = std::ceil(imageCoordinates.x);

    
    float y_low = std::floor(imageCoordinates.y);
    float y_high = std::ceil(imageCoordinates.y);


    glm::vec3 color_low_left = m_pixels[linearize2DCoordinates((unsigned int) x_low, (unsigned int) y_low)];
    glm::vec3 color_low_right = m_pixels[linearize2DCoordinates((unsigned int) x_high, (unsigned int) y_low)];
    glm::vec3 color_high_left = m_pixels[linearize2DCoordinates((unsigned int) x_low, (unsigned int) y_high)];
    glm::vec3 color_high_right = m_pixels[linearize2DCoordinates((unsigned int) x_high, (unsigned int) y_high)];

    glm::vec3 color_low = linearInterpolation(x_low, x_high, color_low_left, color_low_right, imageCoordinates.x);
    glm::vec3 color_high = linearInterpolation(x_low, x_high, color_high_left, color_high_right, imageCoordinates.x);

    glm::vec3 color = linearInterpolation(y_low, y_high, color_low, color_high, imageCoordinates.y);

    return color;
}

// interpolates linearly between 2 colors
glm::vec3 Image::linearInterpolation(float low, float high, glm::vec3 color_low, glm::vec3 color_high, float p) const {
    float c = (p - low) / (high - low);

    return (1 - c) * color_low + c * color_high;
}

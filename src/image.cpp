#include "image.h"
#include "disable_all_warnings.h"
DISABLE_WARNINGS_PUSH()
#include <stb_image.h>
#include <stb_image_write.h>
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

    for (size_t i = 0; i < m_width * m_height * numChannels; i += numChannels) {
        m_pixels.emplace_back(pixels[i + 0] / 255.0f, pixels[i + 1] / 255.0f, pixels[i + 2] / 255.0f);
    }

    stbi_image_free(pixels);

    // add the current texture to the mipmap vector
    _mipmap.push_back(m_pixels);
    curr_height = m_height;
    curr_width = m_width;


    //if the mipmap can be built, build it
    if(canUseMipmapping())
        initMipmap();

    //used for debugging
    //printMipmap();
}

glm::vec3 Image::getPixel(const glm::vec2& textureCoordinates, float pixelDensity) const
{
    // deal with the out of bounds cases

    //std::cout << pixelDensity << std::endl;

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
        return nearestNeighbor(imageCoordinates, 0, m_width, m_height);
    else if(_filteringMethod == TextureFiltering::Bilinear)
        return bilinearInterpolation(imageCoordinates, 0, m_width, m_height);
    else if(_filteringMethod == TextureFiltering::MipMappingNearestLevelNearestNeighbor)
        return nearestLevelMipmapping(inBoundsTextureCoordinates, pixelDensity);
    else if(_filteringMethod == TextureFiltering::MipMappingNearestLevelBilinear)
        return nearestLevelBilinear(inBoundsTextureCoordinates, pixelDensity);
    else if(_filteringMethod == TextureFiltering::Trilinear)
        return trilinearInterpolation(inBoundsTextureCoordinates, pixelDensity);

    return glm::vec3(1);
    
}

// the texture coordinates are in the space [0, 1] x [0, 1], with (0, 0) being the lower-left corner, and (1, 1) the upper right corner
// this functions maps them to the space [0, m_width - 1] x [0, m_height - 1]
glm::vec2 Image::toImageCoordinates(glm::vec2 textureCoordinates, unsigned int level) const {

    if(!levelIsValid(level)) {
        std::cerr << "toImageCoordinates: Invalid level";
        return glm::vec2(0);
    }

    unsigned int w, h;
    if(!getWidthHeightForLevel(w, h, level)) {
        std::cerr << "toImageCoordinates: Could not get width and height for the specified level" << std::endl;
    }

    // the y value must be inverted, since the origin is at the top of the screen for an image
    return glm::vec2(textureCoordinates.x * (w - 1), (1.0f - textureCoordinates.y) * (h - 1));
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
unsigned int Image::linearize2DCoordinates(unsigned int x, unsigned int y, unsigned int width) const {
    // x is the column and y is the row!
    return y * width + x;
}   

// returns the corresponding texture color after performing the nearest neighbor approximation
glm::vec3 Image::nearestNeighbor(glm::vec2 imageCoordinates, int level, unsigned int width, unsigned int height) const
{
    unsigned int x = round(imageCoordinates.x);
    unsigned int y = round(imageCoordinates.y);

    // correction, in case some error has occured, to avoid index out of bounds
    if(x >= width) {
        x = width - 1;
        std::cerr << "Nearest Neighbor: The obtained x was larger than the image width" << std::endl;
    }

    if (x < 0) {
        x = 0;
        std::cerr << "Nearest Neighbor: The obtained x was smaller than 0" << std::endl;
    }

    if (y >= height) {
        y = height - 1;
        std::cerr << "Nearest Neighbor: The obtained y was larger than the image height" << std::endl;
    }

    if (y < 0) {
        y = 0;
        std::cerr << "Nearest Neighbor: The obtained y was smaller than 0" << std::endl;
    }

    return _mipmap[level][linearize2DCoordinates(x, y, width)];
}

// returns the corresponding texture color after applying bilinear interpolation
glm::vec3 Image::bilinearInterpolation(glm::vec2 imageCoordinates, int level, unsigned int width, unsigned int height) const
{
    float x_low = std::floor(imageCoordinates.x);
    float x_high = std::ceil(imageCoordinates.x);

    
    float y_low = std::floor(imageCoordinates.y);
    float y_high = std::ceil(imageCoordinates.y);


    glm::vec3 color_low_left = _mipmap[level][linearize2DCoordinates((unsigned int) x_low, (unsigned int) y_low, width)];
    glm::vec3 color_low_right = _mipmap[level][linearize2DCoordinates((unsigned int)x_high, (unsigned int)y_low, width)];
    glm::vec3 color_high_left = _mipmap[level][linearize2DCoordinates((unsigned int)x_low, (unsigned int)y_high, width)];
    glm::vec3 color_high_right = _mipmap[level][linearize2DCoordinates((unsigned int)x_high, (unsigned int)y_high, width)];

    glm::vec3 color_low = linearInterpolation(x_low, x_high, color_low_left, color_low_right, imageCoordinates.x);
    glm::vec3 color_high = linearInterpolation(x_low, x_high, color_high_left, color_high_right, imageCoordinates.x);

    glm::vec3 color = linearInterpolation(y_low, y_high, color_low, color_high, imageCoordinates.y);

    return color;
}

// makes use of nearest level mipmapping to get the color based on the image coordinates.
// in case of an error, such as the mipmap not being initialized, returns a white pixel
glm::vec3 Image::nearestLevelMipmapping(glm::vec2 textureCoordinates, float pixelDensity) const {
    if(!isMipmapInit())
        return glm::vec3(1);


    unsigned int best_level;
    bool succ = getBestLevelMipmap(best_level, pixelDensity, 0);

    // returns the image coordinates for this level
    glm::vec2 imageCoordinates = toImageCoordinates(textureCoordinates, best_level);

    if(!succ)
        return glm::vec3(1);

    unsigned int w_best, h_best;
    succ = getWidthHeightForLevel(w_best, h_best, best_level);

    if(!succ)
        return glm::vec3(1);

    return nearestNeighbor(imageCoordinates, best_level, w_best, h_best);
}

// same as nearestLevelMipmapping but with bilinear interpolation instead
glm::vec3 Image::nearestLevelBilinear(glm::vec2 textureCoordinates, float pixelDensity) const
{
    if (!isMipmapInit())
        return glm::vec3(1);

    unsigned int best_level;
    bool succ = getBestLevelMipmap(best_level, pixelDensity, 0);

    // returns the image coordinates for this level
    glm::vec2 imageCoordinates = toImageCoordinates(textureCoordinates, best_level);

    if (!succ)
        return glm::vec3(1);

    unsigned int w_best, h_best;
    succ = getWidthHeightForLevel(w_best, h_best, best_level);

    if (!succ)
        return glm::vec3(1);

    return bilinearInterpolation(imageCoordinates, best_level, w_best, h_best);
}

// performs mipmapping with trilinear interpolation
glm::vec3 Image::trilinearInterpolation(glm::vec2 textureCoordinates, float pixelDensity) const
{
    if (!isMipmapInit())
        return glm::vec3(1);

    unsigned int best_level_lower;
    unsigned int best_level_higher;

    bool succ1 = getBestLevelMipmap(best_level_higher, pixelDensity, 1);
    bool succ2 = getBestLevelMipmap(best_level_lower, pixelDensity, 2);

    if(!succ1 && !succ2)
        return glm::vec3(1);

    // if we don't have to mipmaps to interpolate our values in between, then it's just as nearestLevelBilinear
    if(!succ1) {
        glm::vec2 imageCoordinates = toImageCoordinates(textureCoordinates, best_level_higher);
        return nearestLevelBilinear(imageCoordinates, pixelDensity);
    }

    if (!succ2)
    {
        glm::vec2 imageCoordinates = toImageCoordinates(textureCoordinates, best_level_lower);
        return nearestLevelBilinear(imageCoordinates, pixelDensity);
    }

   

    unsigned int w_low, h_low, w_high, h_high;
    if (!getWidthHeightForLevel(w_low, h_low, best_level_lower))
        return glm::vec3(1);


    if (!getWidthHeightForLevel(w_high, h_high, best_level_higher))
        return glm::vec3(1);

    float ratio_low = computePixelRatio(best_level_lower, pixelDensity);
    float ratio_high = computePixelRatio(best_level_higher, pixelDensity);

    glm::vec2 imageCoordinates_lower = toImageCoordinates(textureCoordinates, best_level_lower);
    glm::vec2 imageCoordinates_higher = toImageCoordinates(textureCoordinates, best_level_higher);

    glm::vec3 color_low = bilinearInterpolation(imageCoordinates_lower, best_level_lower, w_low, h_low);
    glm::vec3 color_high = bilinearInterpolation(imageCoordinates_higher, best_level_higher, w_high, h_high);

    // linear interpolate between the low and the high result
    glm::vec3 color = linearInterpolation(ratio_low, ratio_high, color_low, color_high, 1.0f);

    return color;
}

// interpolates linearly between 2 colors
glm::vec3 Image::linearInterpolation(float low, float high, glm::vec3 color_low, glm::vec3 color_high, float p) const {
    float c = (p - low) / (high - low);

    return (1 - c) * color_low + c * color_high;
}


// METHODS USED IN MIPMAPPING

// creates a new texture of reduced resolution corresponding to the next level in the mipmap hierarchy
void Image::getReducedResolutionTexture(unsigned int or_width, unsigned int or_height, std::vector<glm::vec3> &original, unsigned int &re_width, unsigned int &re_height, std::vector<glm::vec3> &reduced) const
{

    re_height = or_height / 2;
    re_width = or_width / 2;

    
    for (int x = 0, re_x = 0; x + 1 < or_width && re_x < re_width; x += 2, re_x++) {
        for (int y = 0, re_y = 0; y + 1 < or_height && re_y < re_height; y += 2, re_y++) {
            int left_upper = linearize2DCoordinates(x, y, or_width);
            int left_lower = linearize2DCoordinates(x, y + 1, or_width);
            int right_upper = linearize2DCoordinates(x + 1, y, or_width);
            int right_lower = linearize2DCoordinates(x + 1, y + 1, or_width);

            glm::vec3 left_upper_color = original[left_upper];
            glm::vec3 left_lower_color = original[left_lower];
            glm::vec3 right_upper_color = original[right_upper];
            glm::vec3 right_lower_color = original[right_lower];


            reduced[linearize2DCoordinates(re_x, re_y, re_width)] = 0.25f * (left_upper_color  + left_lower_color + right_upper_color + right_lower_color);
        }
    }
}

// Returns true if the provided texture supports mipmapping, i.e. if the texture is 2^N x 2^M
bool Image::canUseMipmapping() const {
    return ((m_height & (m_height - 1)) == 0) && ((m_width & (m_width - 1)) == 0) && (m_width == m_height);
}

// initializes the mipmap if it is possible and if it wasn't already initialized
void Image::initMipmap() {
    if(!canUseMipmapping())
        return;

    if(isMipmapInit())
        return;

    _mipmap_init = true;



    int k = m_width;
    unsigned int w = m_width, h = m_height;

    while(k > 1) {

        std::vector<glm::vec3> current(w * h);

        getReducedResolutionTexture(k, k, _mipmap[_mipmap.size() - 1], w, h, current);
        _mipmap.push_back(current);

        k /= 2;
       
    }

}

// Returns true if the mipmap is initialized, otherwise false
bool Image::isMipmapInit() const {
    return _mipmap_init;
}

// returns the number of mipmap levels,
int Image::getNumMipmapLevels() const {
    return _mipmap.size();
}

// switches to another mipmap level. Returns false and does not do anything if it was unsuccesful, i.e. if the mipmap was not initialized or the level is invalid
bool Image::switchToLevel(unsigned int level) {
    if(!isMipmapInit())
        return false;


    if(!levelIsValid(level))
        return false;


    bool succ = getWidthHeightForLevel(curr_width, curr_height, level);

    if(!succ)
        return false;

    curr_level = level;

    curr_texture = _mipmap[level];

    return true;
}

// returns true if the level is valid, otherwise false.
bool Image::levelIsValid(unsigned int level) const {
    return (level >= 0) && (level < getNumMipmapLevels());
}

// returns the currently used level of the mipmap. 
int Image::getCurrMipmapLevel() const {
    return curr_level;
}

// returns the width and the height corresponding to the provided level. If the level is invalid,
// the width and the height remain unchanged
bool Image::getWidthHeightForLevel(unsigned int &width, unsigned int &height, unsigned int level) const{
    if(!levelIsValid(level))
        return false;

    width = m_width >> level;
    height = m_height >> level;

    return true;
}

// compute the ratio between the number of screen pixels and the number of texture pixels for the provided mipmap level
// returns a negative value in the case of an invalid level or uninitialized mipmap, or if the operation is unsuccesful for some other reasons.
float Image::computePixelRatio(unsigned int level, float pixelDensity) const {
    if(!isMipmapInit())
        return -1.0f;

    if(!levelIsValid(level))
        return -1.0f;

    unsigned int w, h;
    bool succ = getWidthHeightForLevel(w, h, level);

    if(!succ)
        return -1.0f;

    if(w * h == 0)
        return -1.0f;

    return (pixelDensity) / (float)(w * h);
}



// returns the best mipmap level with a higher resolution than the screen resolution.
bool Image::getBestLevelMipmap(unsigned int &best_level, float pixelDensity, unsigned int mode = 0) const {
    if(!isMipmapInit()) {
        std::cerr << "Warning! Mipmap not initialized!" << std::endl;
        return false;
    }
    
    std::vector<std::pair<float, unsigned int>> distance_from_1;

    for(int level = 0; level < getNumMipmapLevels(); level++) {
        float ratio = computePixelRatio(level, pixelDensity);
        if(mode == 0 || (mode == 1 && ratio >= 1) || (mode == 2 && ratio <= 1)) {
            float dist = std::abs(1 - ratio);
            distance_from_1.push_back(std::make_pair(dist, level));
        }
    }

    if(distance_from_1.size() == 0)
        return false;

    std::sort(distance_from_1.begin(), distance_from_1.end());

    best_level = distance_from_1[0].second;

    return true;
}


// used only for debugging purposes
void Image::printMipmap() {

    

    for(int level = 0; level < getNumMipmapLevels(); level++) {
       switchToLevel(level);
        std::string filename = "level_";
        filename.push_back((char)(level + '0'));
        filename += ".bmp";

        std::vector<glm::u8vec4> textureData8Bits(curr_texture.size());
        std::transform(std::begin(curr_texture), std::end(curr_texture), std::begin(textureData8Bits),
                       [](const glm::vec3 &color) {
                           const glm::vec3 clampedColor = glm::clamp(color, 0.0f, 1.0f);
                           return glm::u8vec4(glm::vec4(clampedColor, 1.0f) * 255.0f);
                       });


        stbi_write_bmp(filename.c_str(), (int) curr_width, (int) curr_height, 4, textureData8Bits.data());
    }
}
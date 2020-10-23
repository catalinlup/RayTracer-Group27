#pragma once
#include "disable_all_warnings.h"
#include <vector>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <iostream>

enum class FilteringOption
{
    None,
    Bloom,
    BloomWithReinhardHdr,
    BloomWithExposureHdr,
    OnlyLight,
    OnlyLightWithKernel

};

enum class Kernel {
    BoxKernel
};

class Screen {
public:
    Screen(const glm::ivec2& resolution);

    void clear(const glm::vec3& color);
    void setPixel(int x, int y, const glm::vec3& color);

    void writeBitmapToFile(const std::filesystem::path& filePath);



    void draw();



    // controls weather or not the bloom filter should be used on the live output or not
    void setBloomFilterLive(bool bloomFilterLive);

    // Used in bloom filtering
    // set to true if the bloom filter should be active, false otherwise.
    // the default no bloom filtering
    void setBloomFilter(FilteringOption option);

    // sets the kernel to be used on the light image in the process of bloom filtering
    // the default kernel is a box kernel
    void setKernel(Kernel kernel);

    // configures the gamma value used in the gamma correction
    void setGammaValue(float gamma);

    // set to true to enable gamma correction, false to disable it
    void enableGammaCorrection(float gammaCorrection);

    // set the exposure level for the exposure bloom filtering
    // the default is 0.5
    void setExposure(float exposure);

    // set the filter size, the default value is 5
    void setFilterSize(int filterSize);


private:
    glm::ivec2 m_resolution;
    std::vector<glm::vec3> m_textureData;

    uint32_t m_texture; 

    bool _bloomFilterLive = false; // controls weather of not the bloom filter should be used in the live preview as well, or only in the final output

    // draws an image of the same size as the screen resolution.
    // used for debugging
    void drawImage(std::vector<glm::vec3> &image);

    // Variables used in bloom filtering

    FilteringOption _filtering_option = FilteringOption::None; // controls the bloom filter, the default is no bloom filtering
    Kernel _kernel = Kernel::BoxKernel; // the default is a box kernel

    float _gamma = 2.2f;

    // true if gamma correction is used, false otherwise
    bool _gammaCorrection = false;

    float _exposure = 0.5;

    int _filterSize = 5;

    // methods used in bloom filtering

    // applies the bloom effect on the texture data, based on the configured settings
    void applyBloomEffect();


    // based on the current texture, returns another texture containing only the pixels with values that are larger than 1
    void filterLightPixels(std::vector<glm::vec3> &lightPixels);

    // applies a box kernel on the provided image, of the provided size
    // the box kernel takes the average of all values that fall into it
    // if the operation failed returns false, otherwise true
    bool applyKernel(Kernel kernel, std::vector<glm::vec3> &image, int width, int height);

    // implements a simple box kernel
    glm::vec3 boxKernel(int x, int y, int filterSize, std::vector<glm::vec3> &image, int width, int height);

    // add the values of the 2 images
    // if the operation failed returns false, otherwise true
    bool addImages(std::vector<glm::vec3>& image1, std::vector<glm::vec3>& image2, std::vector<glm::vec3>& addedImage);

    // clamp the values of the pixels that exceed 1
    glm::vec3 clamp(glm::vec3 pixel);

    // convert HDR values to Non-HDR values using Reinhard's algorithm
    glm::vec3 reinhardToneMap(glm::vec3 pixel);

    // convert HDR values to Non-HDR values using Exposure Tone Mapping
    glm::vec3 exposureToneMap(glm::vec3 pixel);

    // applies gamma correction to the final pixels
    glm::vec3 gammaCorrection(glm::vec3 pixel);

    // convert pixel to grayscale in order to get the 
    float convertToGrayscale(glm::vec3 pixel);

    // returns the pixel at the corresponding coordinates
    // if the pixels are out of bounds, returns a black pixel (border)
    glm::vec3 getPixel(int x, int y, std::vector<glm::vec3>& image);


};

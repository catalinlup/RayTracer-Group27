#include "screen.h"
#include "disable_all_warnings.h"
#include "opengl_includes.h"
DISABLE_WARNINGS_PUSH()
#include <algorithm>
#include <glm/common.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <stb_image_write.h>
#include <string>
DISABLE_WARNINGS_POP()

Screen::Screen(const glm::ivec2& resolution)
    : m_resolution(resolution)
    , m_textureData(size_t(resolution.x * resolution.y), glm::vec3(0.0f))
{
    // Generate texture
    glGenTextures(1, &m_texture);
    glBindTexture(GL_TEXTURE_2D, m_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void Screen::clear(const glm::vec3& color)
{
    std::fill(std::begin(m_textureData), std::end(m_textureData), color);
}

void Screen::setPixel(int x, int y, const glm::vec3& color)
{
    // In the window/camera class we use (0, 0) at the bottom left corner of the screen (as used by GLFW).
    // OpenGL / stbi like the origin / (-1,-1) to be at the TOP left corner so transform the y coordinate.
    const int i = (m_resolution.y - 1 - y) * m_resolution.x + x;
    m_textureData[i] = glm::vec4(color, 1.0f);
}

void Screen::writeBitmapToFile(const std::filesystem::path& filePath)
{

    std::vector<glm::u8vec4> textureData8Bits(m_textureData.size());
    std::transform(std::begin(m_textureData), std::end(m_textureData), std::begin(textureData8Bits),
        [](const glm::vec3& color) {
            const glm::vec3 clampedColor = glm::clamp(color, 0.0f, 1.0f);
            return glm::u8vec4(glm::vec4(clampedColor, 1.0f) * 255.0f);
        });

    std::string filePathString = filePath.string();
    stbi_write_bmp(filePathString.c_str(), m_resolution.x, m_resolution.y, 4, textureData8Bits.data());
}

// postProcesses, i.e. applies gamma correction or bloom effects
void Screen::postprocessImage() {
    if (_bloomFilterLive)
    {
        applyBloomEffect();
    }

    // use gamma correction if enabled
    // if gamma correction is set, apply gamma correction
    if (_gammaCorrection)
    {
        for (int i = 0; i < m_textureData.size(); i++)
            m_textureData[i] = gammaCorrection(m_textureData[i]);
    }
}

void Screen::draw()
{

    drawImage(m_textureData);

    // glPushAttrib(GL_ALL_ATTRIB_BITS);

    // glBindTexture(GL_TEXTURE_2D, m_texture);
    // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, m_resolution.x, m_resolution.y, 0, GL_RGB, GL_FLOAT, m_textureData.data());

    // glDisable(GL_LIGHTING);
    // glDisable(GL_LIGHT0);
    // glDisable(GL_COLOR_MATERIAL);
    // glDisable(GL_NORMALIZE);
    // glColor3f(1.0f, 1.0f, 1.0f);

    // glEnable(GL_TEXTURE_2D);
    // glActiveTexture(GL_TEXTURE0);
    // glBindTexture(GL_TEXTURE_2D, m_texture);

    // glMatrixMode(GL_MODELVIEW);
    // glPushMatrix();
    // glLoadIdentity();
    // glMatrixMode(GL_PROJECTION);
    // glPushMatrix();
    // glLoadIdentity();

    // glBegin(GL_QUADS);
    // glTexCoord2f(0.0f, 1.0f);
    // glVertex3f(-1.0f, -1.0f, 0.0f);
    // glTexCoord2f(1.0f, 1.0f);
    // glVertex3f(+1.0f, -1.0f, 0.0f);
    // glTexCoord2f(1.0f, 0.0f);
    // glVertex3f(+1.0f, +1.0f, 0.0f);
    // glTexCoord2f(0.0f, 0.0f);
    // glVertex3f(-1.0f, +1.0f, 0.0f);
    // glEnd();

    // glMatrixMode(GL_PROJECTION);
    // glPopMatrix();
    // glMatrixMode(GL_MODELVIEW);
    // glPopMatrix();

    // glPopAttrib();
}

// draws an image of the same size as the screen resolution.
// used for debugging
void Screen::drawImage(std::vector<glm::vec3> &image) {

    if(m_resolution.x * m_resolution.y != image.size()) {
        std::cerr << "Could not draw image! The provided image does not match the screen resolution!" << std::endl;
        return;
    }

   


    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glBindTexture(GL_TEXTURE_2D, m_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, m_resolution.x, m_resolution.y, 0, GL_RGB, GL_FLOAT, image.data());

    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_NORMALIZE);
    glColor3f(1.0f, 1.0f, 1.0f);

    glEnable(GL_TEXTURE_2D);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, m_texture);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 1.0f);
    glVertex3f(-1.0f, -1.0f, 0.0f);
    glTexCoord2f(1.0f, 1.0f);
    glVertex3f(+1.0f, -1.0f, 0.0f);
    glTexCoord2f(1.0f, 0.0f);
    glVertex3f(+1.0f, +1.0f, 0.0f);
    glTexCoord2f(0.0f, 0.0f);
    glVertex3f(-1.0f, +1.0f, 0.0f);
    glEnd();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glPopAttrib();
}

// methods used in bloom filtering

// controls weather or not the bloom filter should be used on the live output or not
void Screen::setBloomFilterLive(bool bloomFilterLive) {
    _bloomFilterLive = bloomFilterLive;
}

// Used in bloom filtering
// set to true if the bloom filter should be active, false otherwise.
// the default no bloom filtering
void Screen::setBloomFilter(FilteringOption option) {
    _filtering_option = option;
}

// sets the kernel to be used on the light image in the process of bloom filtering
// the default kernel is a box kernel
void Screen::setKernel(Kernel kernel) {
    _kernel = kernel;
}

// sets the number of times the kernel should be applied. The default is 1
void Screen::setKernelNumRepetitions(int repetitions) {
    _kernel_num_repetitions = repetitions;
    _kernel_num_repetitions = glm::max(1, _kernel_num_repetitions);
}

// configures the gamma value used in the gamma correction
void Screen::setGammaValue(float gamma) {
    _gamma = gamma;
}

// set to true to enable gamma correction, false to disable it
void Screen::enableGammaCorrection(float gammaCorrection) {
    _gammaCorrection = gammaCorrection;
}

// set the exposure level for the exposure bloom filtering
// the default is 0.5
void Screen::setExposure(float exposure) {
    _exposure = exposure;
}

// sets the sigma value to be used in the gaussian blurr (default 2.0)
void Screen::setSigma(float sigma) {
    _sigma = sigma;
    _sigma = glm::max(0.001f, _sigma);
}

// set the box filter size, the default value is 5
void Screen::setFilterSize(int filterSize) {
    _filterSize = filterSize;
    filterSize = glm::max(1, filterSize);
}

// applies the bloom effect on the texture data, based on the configured settings
void Screen::applyBloomEffect() {

    if(_filtering_option == FilteringOption::None)
        return;

    if(_filtering_option == FilteringOption::OnlyLight) {
        std::vector<glm::vec3> lightImage;
        filterLightPixels(lightImage);
        m_textureData.clear();
        m_textureData = lightImage;
        return;
    }

    if(_filtering_option == FilteringOption::OnlyLightWithKernel) {
        std::vector<glm::vec3> lightImage;
        filterLightPixels(lightImage);
        applyKernel(_kernel, lightImage, m_resolution.x, m_resolution.y);
        m_textureData.clear();
        m_textureData = lightImage;
        return;
    }

    // can the parts of the image with brightness that is larger than 1
    std::vector<glm::vec3> lightImage;
    filterLightPixels(lightImage);

    // apply the kernel the configured number of times
    for(int i = 1; i <= _kernel_num_repetitions; i++)
        applyKernel(_kernel, lightImage, m_resolution.x, m_resolution.y);

    std::vector<glm::vec3> data_copied(m_textureData);
    addImages(data_copied, lightImage, m_textureData);

    // apply tone mapping and gamma correction
    for(int i = 0; i < m_textureData.size(); i++) {
        if(_filtering_option == FilteringOption::Bloom)
            m_textureData[i] = clamp(m_textureData[i]);
        else if(_filtering_option == FilteringOption::BloomWithReinhardHdr)
            m_textureData[i] = reinhardToneMap(m_textureData[i]);
        else if(_filtering_option == FilteringOption::BloomWithExposureHdr)
            m_textureData[i] = exposureToneMap(m_textureData[i]);

  
    }
}


// based on the current texture, returns another texture containing only the pixels with values that are larger than 1
void Screen::filterLightPixels(std::vector<glm::vec3> &lightPixels) {

    lightPixels.clear();

    for(auto pixel : m_textureData) {
        float brightness = convertToGrayscale(pixel);

        if(brightness >= 1.0f)
            lightPixels.push_back(pixel);
        else
            lightPixels.push_back(glm::vec3(0));
    }
}

// applies a box kernel on the provided image, of the provided size
// the box kernel takes the average of all values that fall into it
bool Screen::applyKernel(Kernel kernel, std::vector<glm::vec3> &image, int width, int height) {

    if(width * height != image.size()) {
        std::cerr << "applyKernel: Size mismatch" << std::endl;

        return false;
    }

    std::vector<glm::vec3> source(image);

    for(int x = 0; x < width; x++) {
        for(int y = 0; y < height; y++) {

            if (kernel == Kernel::BoxKernel)
                image[y * width + x] = boxKernel(x, y, _filterSize, source, width, height);
            else if (kernel == Kernel::GaussianKernel)
                image[y * width + x] = gaussianKernel(x, y, _filterSize, source, width, height);
        }
    }

    return true;
}

// implements a simple box kernel
glm::vec3 Screen::boxKernel(int x, int y, int filterSize, std::vector<glm::vec3> &image, int width, int height) {

    glm::vec3 sum = glm::vec3(0.0f);

    for(int i = -filterSize; i < filterSize + 1; i++) {
        for(int j = -filterSize; j < filterSize + 1; j++) {
            sum += getPixel(x + i, y + j, image);
        }
    }

    sum /= (2 * filterSize + 1) * (2 * filterSize + 1);

    return sum;
}

// implements a gaussian kernel
glm::vec3 Screen::gaussianKernel(int x, int y, int filterSize, std::vector<glm::vec3> &image, int width, int height) {
    glm::vec3 sum = glm::vec3(0.0f);

    for(int i = -filterSize; i < filterSize + 1; i++) {
        for(int j = -filterSize; j < filterSize + 1; j++) {
            sum += gaussianFunction((float) i, (float) j) * getPixel(x + i, y + j, image);
        }
    }

    return sum;
}

// computes the value of the 2d gaussian function with mean 0 and standard deviation _sigma
float Screen::gaussianFunction(float x, float y) {
    return (1 / (_sigma * _sigma * 2 * M_PI)) * glm::exp(- (x * x + y * y) / (2 * _sigma * _sigma));
}

// add the values of the 2 images of the same size
// if the images don't have the same size, it does not do anything and it returns false
bool Screen::addImages(std::vector<glm::vec3> &image1, std::vector<glm::vec3> &image2, std::vector<glm::vec3> &addedImage) {
    if(image1.size() != image2.size()) {
        std::cerr << "addImages: Size mismatch" << std::endl;
        return false;
    }

    addedImage.clear();

    for(int i = 0; i < image1.size(); i++) {
        addedImage.push_back(image1[i] + image2[i]);
    }
}

// clamp the values of the pixels that exceed 1
glm::vec3 Screen::clamp(glm::vec3 pixel) {
    return glm::clamp(pixel, glm::vec3(0), glm::vec3(1));
}

// convert HDR values to Non-HDR values using Reinhard's algorithm
glm::vec3 Screen::reinhardToneMap(glm::vec3 pixel) {
    return pixel / (pixel + glm::vec3(1.0f));
}

// convert HDR values to Non-HDR values using Exposure Tone Mapping
glm::vec3 Screen::exposureToneMap(glm::vec3 pixel) {
    return glm::vec3(1.0f) - glm::exp(-pixel * _exposure);
}

// applies gamma correction to the final pixels
glm::vec3 Screen::gammaCorrection(glm::vec3 pixel) {
    return glm::pow(pixel, glm::vec3(1.0f / _gamma));
}

// convert pixel to grayscale in order to get the
float Screen::convertToGrayscale(glm::vec3 pixel) {
    return glm::dot(pixel, glm::vec3(0.2126, 0.7152, 0.0722));
}

// returns the pixel at the corresponding coordinates
// if the pixels are out of bounds, returns a black pixel (border)
glm::vec3 Screen::getPixel(int x, int y, std::vector<glm::vec3> &image) {
    if(x < 0 || y < 0 || x >= m_resolution.x || y >= m_resolution.y) {
        return glm::vec3(0);
    }

    return image[y * m_resolution.x + x];
}

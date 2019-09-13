#include "tuple.h"
namespace vgi {

/**
 * Convert a color from <TType> to <TTypeOut>.
 */
template<int TDim, typename TType, typename TTypeOut>
Tuple<TDim, TTypeOut, kColor> convert(const Tuple<TDim, TType, kColor>& source) {
    Tuple<TDim, uint8_t, kColor> out;
    for (int it=0;it<TDim;it++)
        out.data[it] = static_cast<TTypeOut>(source[it]);
    return out;
}

/**
 * Convert an image from <TType> to uint8_t.
 */
template<int TDim, typename TType>
void convert(const Image<TDim, TType>& image, Image<TDim, uint8_t>& out) {
    out.allocate(image.size);
    for (size_t it = 0; it < image.area(); it++) {
        out.data[it] = convert<3, double, uint8_t>(image.data[it] * 255.0);
    }
}

/**
 * Save the image to file using PPM (uncompressed) format.
 */
void save(const Image<3, uint8_t>& image, const char * filename) {
    std::ofstream file;
    file.open(filename);
    file << "P" << 3 << "\n\n";
    file << image.size[0] << " " << image.size[1] << "\n255\n";
    for (size_t it = 0; it < image.area(); it++) {
        file << (int)image.data[it][0] << " " << (int)image.data[it][1] << " " << (int)image.data[it][2] << "\n";
    }
    file.close();
}

} // namespace vgi

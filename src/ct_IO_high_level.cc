#include "ct_IO_high_level.hpp"
#include <zlib.h>

namespace CT
{
  // using zlib to compress searlized data
  // ref: http://hamelot.io/programming/compression-of-simulation-data-using-zlib/
  // we will use GZip from zlib
  
  // read gz file, and convert to string
  void gz_to_string(const std::string in_state_file, std::string & data)
  {
    //open the file for reading in binary mode
    gzFile gz_file = gzopen(in_state_file.c_str(), "rb");

    //this variable will hold the size of the file
    unsigned long int size;
    //we wrote out a unsigned long int when storing the file
    //read this back in to get the size of the uncompressed data
    gzread(gz_file, (void*) &size, sizeof(size));

    //resize the string
    data.resize(size / sizeof(char));
    //read in and uncompress the entire data set at once
    gzread(gz_file, (void*) data.data(), size);
    gzclose(gz_file);
  }

  // convert istringstream to gz file
  void stringstream_to_gz(std::string out_state_file, std::string & data)
  {
    gzFile out_gz_file = gzopen(out_state_file.c_str(), "wb");

    //Write the size of the stream, this is needed so that we know
    //Get the size of the stream

    unsigned long int file_size = data.size();
    //how much to read back in later
    gzwrite(out_gz_file, (void*) &file_size, sizeof(file_size));
    //Write the data
    gzwrite(out_gz_file, (void*) (data.data()), file_size);
    //close the file
    gzclose(out_gz_file);
  }
}

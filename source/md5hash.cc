// *** release/md5hash.cc ***
// Author: Kevin Wolz, date: 08/2018
//
// Generates a hash ID from an arbitrary file whose path is given, using the 
// md5 file hashing algorithm from the OpenSSL library.

#include "md5hash.h"


// Print the MD5 sum as hex-digits.
std::string PrintHashToString(unsigned char* md) {
    std::stringstream ss;
    for(int i=0; i < MD5_DIGEST_LENGTH; i++) {
        ss<<std::to_string(md[i]);
    }
    return ss.str();
}


// Get the size of the file by its file descriptor.
unsigned long GetSizeByFd(int fd) {
    struct stat statbuf;
    if(fstat(fd, &statbuf) < 0) exit(-1);
    return statbuf.st_size;
}


// Outputs the 7-digit string hash ID of the file located at file_path.
std::string HashFile(std::string& file_path) {
    int file_descript;
    unsigned long file_size;
    char* file_buffer;
    unsigned char result[MD5_DIGEST_LENGTH];

    std::cout<<"Creating md5 hash using file: "<<file_path<<std::endl;

    const char* path = file_path.c_str();
    file_descript = open(path, O_RDONLY);
    if(file_descript < 0) exit(-1);

    file_size = GetSizeByFd(file_descript);
    std::cout<<"file size: "<<file_size<<std::endl;

    file_buffer = static_cast<char*>(mmap((caddr_t)0, file_size, PROT_READ,
                                     MAP_SHARED, file_descript, 0));
    MD5((unsigned char*) file_buffer, file_size, result);
    munmap(file_buffer, file_size); 

    std::string hash = PrintHashToString(result);
    std::cout<<"hash = "<<hash<<std::endl;
    hash.resize(7);
    return hash;
}
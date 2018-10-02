#ifndef SOURCE_MD5HASH_H_
#define SOURCE_MD5HASH_H_

// *** source/md5hash.h ***
// Author: Kevin Wolz, date: 08/2018
//
// Generates a hash ID from an arbitrary file whose path is given, using the 
// md5 file hashing algorithm from the OpenSSL library.

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <iostream>
#include <sstream> 
#include <string>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>

#include <openssl/md5.h>

std::string PrintHashToString(unsigned char* md);
unsigned long GetSizeByFd(int fd);
std::string HashFile(std::string& file_path);


#endif // SOURCE_MD5HASH_H_
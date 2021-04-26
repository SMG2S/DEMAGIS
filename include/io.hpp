#include <iostream>
#include <fstream>
#include <complex>
#include <iomanip>
#include <typeinfo>

template <typename T>
void wrtMatIntoBinary(T *H, std::string path_out, std::size_t size){
  std::ostringstream problem(std::ostringstream::ate);
  problem << path_out;

  std::cout << "]> writing matrix into ";
  std::cout << problem.str();
  std::cout << " of size = " << size << std::endl;

  auto outfile = std::fstream(problem.str().c_str(), std::ios::out | std::ios::binary);

  outfile.write((char*)&H[0], size * sizeof(T));

  outfile.close();

}

template <typename T>
void readMatFromBinary(T *H, std::string path_in, std::size_t size){
  std::ostringstream problem(std::ostringstream::ate);
  problem << path_in;

  std::cout << "]> start reading matrix from binary file ";
  std::cout << problem.str();
  std::cout << " of size = " << size << std::endl;

  std::ifstream infile(problem.str().c_str(), std::ios::binary);

  infile.read((char*)H, sizeof(T) * size);

  infile.close();
}



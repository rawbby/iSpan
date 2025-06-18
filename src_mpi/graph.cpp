#include "./graph.h"
#include "./wtime.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <system_error>
#include <vector>

template<class T>
static void
read_binary(const std::filesystem::path& file, std::vector<T>& dst)
{
  std::ifstream in(file, std::ios::binary);
  if (!in)
    throw std::system_error(errno, std::generic_category(), "open " + file.string());
  in.read(reinterpret_cast<char*>(dst.data()), static_cast<std::streamsize>(dst.size() * sizeof(T)));
  if (!in)
    throw std::runtime_error("read " + file.string());
}

graph::graph(const std::filesystem::path& fw_beg_file,
             const std::filesystem::path& fw_csr_file,
             const std::filesystem::path& bw_beg_file,
             const std::filesystem::path& bw_csr_file)
  : src_count(0)
{
  using namespace std::filesystem;
  const double tm = wtime();

  vert_count = file_size(fw_beg_file) / sizeof(index_t) - 1;
  edge_count = file_size(fw_csr_file) / sizeof(vertex_t);

  fw_beg_pos.resize(vert_count + 1);
  fw_csr.resize(edge_count);
  bw_beg_pos.resize(vert_count + 1);
  bw_csr.resize(edge_count);

  read_binary(fw_beg_file, fw_beg_pos);
  read_binary(fw_csr_file, fw_csr);
  read_binary(bw_beg_file, bw_beg_pos);
  read_binary(bw_csr_file, bw_csr);

  std::cout << "Graph load (success): " << vert_count << " verts, "
            << edge_count << " edges in "
            << wtime() - tm << " s\n";
}

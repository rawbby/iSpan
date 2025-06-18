#pragma once

#include "util.h"
#include <vector>

class graph
{
public:
  std::vector<index_t> fw_beg_pos;
  std::vector<vertex_t> fw_csr;
  std::vector<index_t> bw_beg_pos;
  std::vector<vertex_t> bw_csr;
  std::vector<path_t> weight;
  std::vector<vertex_t> src_list;

  index_t src_count;
  index_t vert_count;
  index_t edge_count;

public:
  graph() = default;
  ~graph() = default;
  graph(const char* fw_beg_file,
        const char* fw_csr_file,
        const char* bw_beg_file,
        const char* bw_csr_file);
  void gen_src() {};
  void groupby() {};
};

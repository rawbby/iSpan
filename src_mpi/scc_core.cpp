#include "color_propagation.h"
#include "fw_bw.h"
#include "graph.h"
#include "scc_common.h"
#include "trim_1_gfq.h"
#include "wtime.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#define INF (-1)

void
scc_detection(
  const graph* g,
  int alpha,
  index_t thread_count,
  std::vector<double>& avg_time,
  int world_rank,
  int world_size,
  int run_time,
  std::vector<vertex_t>& assignment)
{
  std::vector<index_t> max_pivot_list(thread_count);
  std::vector<index_t> max_degree_list(thread_count);

  std::vector<index_t> thread_bin(thread_count);

  std::vector<index_t> front_comm(world_size);
  std::vector<long_t> work_comm(world_size);

  std::vector<bool> color_change(thread_count);

  index_t s = g->vert_count / 32; // todo parameterize especially for tests
  if (g->vert_count % 32 != 0)
    s += 1;
  index_t t = s / world_size;
  if (s % world_size != 0)
    t += 1;
  index_t virtual_count = t * world_size * 32;

  index_t step = t * 32;
  auto vert_beg = std::min<index_t>(step * world_rank, g->vert_count);
  auto vert_end = std::min<index_t>(step * (world_rank + 1), g->vert_count);

  std::vector<unsigned int> sa_compress(s);

  std::vector<index_t> small_queue(virtual_count + 1);
  std::vector<index_t> wcc_fq(virtual_count + 1);
  std::vector<vertex_t> vert_map(g->vert_count + 1);
  std::vector<vertex_t> sub_fw_beg(g->vert_count + 1);
  std::vector<vertex_t> sub_fw_csr(g->edge_count + 1);
  std::vector<vertex_t> sub_bw_beg(g->vert_count + 1);
  std::vector<vertex_t> sub_bw_csr(g->edge_count + 1);

  depth_t* fw_sa;
  depth_t* bw_sa;

  if (posix_memalign((void**)&fw_sa, getpagesize(), sizeof(depth_t) * (virtual_count + 1)))
    perror("posix_memalign");

  if (posix_memalign((void**)&bw_sa, getpagesize(), sizeof(depth_t) * (virtual_count + 1)))
    perror("posix_memalign");
  std::vector<vertex_t> scc_id(virtual_count + 1);

  std::vector<vertex_t> fq_comm(virtual_count + 1);
  std::vector<vertex_t> scc_id_mice(virtual_count + 1);
  std::vector<vertex_t> fw_sa_temp(virtual_count + 1);
  std::vector<vertex_t> color(virtual_count + 1);

  for (long_t i = 0; i < virtual_count + 1; ++i) {

    fw_sa[i] = -1;

    bw_sa[i] = -1;
    scc_id[i] = 0;
    scc_id_mice[i] = -1;
  }

  double start_time = wtime();
  double end_time;

  {
    double time_size_1_first;
    double time_size_1;
    double time_fw;
    double time_bw;
    double time_size_2;
    double time_size_3;
    double time_gfq;
    double pivot_time;
    double time_wcc;
    double time_mice_fw_bw;

    double time_comm;

    MPI_Barrier(MPI_COMM_WORLD);
    double time = wtime();

    trim_1_first(scc_id,
                 g->fw_beg_pos,
                 g->bw_beg_pos,
                 vert_beg,
                 vert_end);
    std::cout << world_rank << ",Computing size_1_first cost," << (wtime() - time) * 1000 << " ms\n";

    double temp_time = wtime();
    MPI_Allgather(MPI_IN_PLACE,
                  0,
                  MPI_INT,
                  scc_id.data(),
                  step,
                  MPI_INT,
                  MPI_COMM_WORLD);

    double time_comm_trim_1 = wtime() - temp_time;
    std::cout << world_rank << ",trim-1 comm time," << time_comm_trim_1 * 1000 << ",ms\n";

    if (world_rank == 0) {
      time_size_1_first = wtime() - time;
    }

    time = wtime();
    vertex_t root = pivot_selection(scc_id,
                                    g->fw_beg_pos,
                                    g->bw_beg_pos,
                                    0,
                                    g->vert_count,
                                    g->fw_csr,
                                    g->bw_csr,
                                    max_pivot_list,
                                    max_degree_list,
                                    world_rank,
                                    thread_count);

    pivot_time = wtime() - time;

    time = wtime();
    fw_bfs(scc_id,
           g->fw_beg_pos,
           g->bw_beg_pos,

           vert_beg,
           vert_end,
           g->fw_csr,
           g->bw_csr,
           fw_sa,
           front_comm,
           root,
           world_rank,
           alpha,
           g->edge_count,
           g->vert_count,
           world_size,
           world_rank,
           step,
           fq_comm,
           sa_compress,
           virtual_count);

    time_fw = wtime() - time;
    for (vertex_t i = 0; i < virtual_count / 32; ++i) {
      sa_compress[i] = 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    time = wtime();
    bw_bfs(scc_id,
           g->fw_beg_pos,
           g->bw_beg_pos,
           vert_beg,
           vert_end,
           g->fw_csr,
           g->bw_csr,
           fw_sa,
           bw_sa,
           front_comm,
           work_comm,
           root,
           world_rank,
           alpha,
           g->edge_count,
           g->vert_count,
           world_size,
           world_rank,
           step,
           fq_comm,
           sa_compress);

    time_bw = wtime() - time;

    printf("bw_bfs done\n");

    time = wtime();

    trim_1_normal(scc_id,
                  g->fw_beg_pos,
                  g->bw_beg_pos,
                  vert_beg,
                  vert_end,
                  g->fw_csr,
                  g->bw_csr);

    if (world_rank == 0) {
      time_size_1 += wtime() - time;
    }

    time = wtime();
    trim_1_normal(scc_id,
                  g->fw_beg_pos,
                  g->bw_beg_pos,
                  vert_beg,
                  vert_end,
                  g->fw_csr,
                  g->bw_csr);

    time_size_1 += wtime() - time;

    time = wtime();

    MPI_Allreduce(MPI_IN_PLACE,
                  scc_id.data(),
                  scc_id.size(),
                  MPI_LONG,
                  MPI_MAX,
                  MPI_COMM_WORLD);

    temp_time = wtime();

    gfq_origin(g->vert_count,
               scc_id,
               small_queue,
               g->fw_beg_pos,
               g->fw_csr,
               g->bw_beg_pos,
               g->bw_csr,
               sub_fw_beg,
               sub_fw_csr,
               sub_bw_beg,
               sub_bw_csr,
               front_comm,
               work_comm,
               world_rank,
               vert_map);

    vertex_t sub_v_count = front_comm[world_rank];

    if (sub_v_count > 0) {
      vertex_t wcc_fq_size = 0;

      for (index_t i = 0; i < sub_v_count; ++i) {
        color[i] = i;
      }

      time = wtime();
      coloring_wcc(
        color,
        sub_fw_beg,
        sub_fw_csr,
        sub_bw_beg,
        sub_bw_csr,
        0,
        sub_v_count);

      if (world_rank == 0) {
        time_wcc += wtime() - time;
      }
      time = wtime();

      process_wcc(0,
                  sub_v_count,
                  wcc_fq,
                  color,
                  wcc_fq_size);

      if (world_rank == 0) {
        printf("color time (ms), %lf, wcc_fq, %lu, time (ms), %lf\n", time_wcc * 1000, wcc_fq_size, 1000 * (wtime() - time));
      }

      mice_fw_bw(color,
                 scc_id_mice,
                 sub_fw_beg,
                 sub_bw_beg,
                 sub_fw_csr,
                 sub_bw_csr,
                 fw_sa_temp,
                 world_rank,
                 world_size,
                 sub_v_count,
                 wcc_fq,
                 wcc_fq_size);

      temp_time = wtime();
      MPI_Allreduce(MPI_IN_PLACE,
                    scc_id_mice.data(),
                    g->vert_count,
                    MPI_LONG,
                    MPI_MAX,
                    MPI_COMM_WORLD);
      time_comm = wtime() - temp_time;

      printf("%d,final comm time,%.3lf\n", world_rank, time_comm * 1000);

      for (index_t i = 0; i < sub_v_count; ++i) {
        vertex_t actual_v = small_queue[i];
        scc_id[actual_v] = small_queue[scc_id_mice[i]];
      }
    }
    if (world_rank == 0) {
      time_mice_fw_bw = wtime() - time;
    }

    if (world_rank == 0 && run_time != 1) {
      avg_time[0] += time_size_1_first + time_size_1 + time_size_2 + time_size_3;

      avg_time[1] += time_fw + time_bw;

      avg_time[2] += time_wcc + time_mice_fw_bw;
      avg_time[4] += time_size_1_first + time_size_1;

      avg_time[5] += time_size_2;
      avg_time[6] += pivot_time;
      avg_time[7] += time_fw;
      avg_time[8] += time_bw;
      avg_time[9] += time_wcc;
      avg_time[10] += time_mice_fw_bw;

      avg_time[13] += time_size_3;
      avg_time[14] += time_gfq;
    }
    if (OUTPUT_TIME) {
      if (world_rank == 0) {
        printf("\ntime size_1_first, %.3lf\ntime size_1, %.3lf\ntime pivot, %.3lf\nlargest fw, %.3lf\nlargest bw, %.3lf\nlargest fw/bw, %.3lf\ntrim size_2, %.3lf\ntrim size_3, %.3lf\nwcc time, %.3lf\nmice fw-bw time, %.3lf\nmice scc time, %.3lf\ntotal time, %.3lf\n", time_size_1_first * 1000, time_size_1 * 1000, pivot_time * 1000, time_fw * 1000, time_bw * 1000, (pivot_time + time_fw + time_bw) * 1000, time_size_2 * 1000, time_size_3 * 1000, time_wcc * 1000, time_mice_fw_bw * 1000, (time_wcc + time_mice_fw_bw) * 1000, (time_size_1_first + time_size_1 + pivot_time + time_fw + time_bw + time_size_2 + time_size_3 + time_wcc + time_mice_fw_bw) * 1000);
      }
    }
#pragma omp barrier
  }
  end_time = wtime() - start_time;
  avg_time[3] += end_time;

  for (auto i = 0; i < g->vert_count; ++i) {
    std::cout << i << ": " << scc_id[i] << "(" << world_rank << ")" << std::endl;
  }
  get_scc_result(scc_id, g->vert_count);

  if (!assignment.empty()) {
    std::unordered_map<vertex_t, std::vector<vertex_t>> scc_components;
    for (vertex_t i = 0; i < g->vert_count; ++i) {
      const auto repr = scc_id[i] == -1 ? i : scc_id[i];
      scc_components[repr].push_back(i);
    }

    for (const auto& [_, vs] : scc_components) {
      const auto min_v = *std::ranges::min_element(vs);
      for (vertex_t v : vs) {
        std::cout << v << ' ';
        assignment[v] = min_v;
      }
      std::cout << std::endl;
    }
  }
}

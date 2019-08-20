#ifndef _UTILS_H_
#define _UTILS_H_

#include <mpi.h>
#include <vector>
#include <string>

void gather_protobuf(
    int local_num, 
    std::vector<std::string> &local_entities, 
    int &total_num,
    std::vector<std::string> &all_entities, 
    int center_node, 
    MPI_Comm comm
  );

void all2all_protobuf(
    std::vector<std::string> &send_entities, 
    std::vector<std::string> &recv_entities,
    MPI_Comm comm
  );

#endif	


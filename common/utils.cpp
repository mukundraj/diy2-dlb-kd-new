#include "utils.h"
#include <cmath>
#include <climits>
#include <algorithm>

using namespace std; 

void gather_protobuf(
    int local_num, 
    vector<string> &local_entities, 
    int &total_num, 
    vector<string> &all_entities,
    int center_node, 
    MPI_Comm comm
  ) 
{
  int comm_size;
  int comm_rank;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  MPI_Allreduce(&local_num, &total_num, 1, MPI_INT, MPI_SUM, comm);

  vector<int> local_length;
  string local_concat;
  for (int i = 0; i < local_num; ++i) {
    string buf = local_entities[i];
    local_length.push_back(buf.size());
    local_concat += buf;
  }

  vector<int> all_length(total_num);
  {
    vector<int> each_num(comm_size);
    if (center_node < 0)
      MPI_Allgather(&local_num, 1, MPI_INT, each_num.data(), 1, MPI_INT, comm);
    else
      MPI_Gather(&local_num, 1, MPI_INT, each_num.data(), 1, MPI_INT, center_node, comm);
    vector<int> displs(comm_size);
    displs[0] = 0;
    for (int i = 1; i < displs.size(); ++i) displs[i] = displs[i - 1] + each_num[i - 1];

    if (center_node < 0)
      MPI_Allgatherv(
          (void*)local_length.data(), local_num, MPI_INT,
	  (void*)all_length.data(), each_num.data(), displs.data(), MPI_INT,
	  comm
        );
    else
      MPI_Gatherv(
          (void*)local_length.data(), local_num, MPI_INT,
	  (void*)all_length.data(), each_num.data(), displs.data(), MPI_INT,
	  center_node, comm
        );
  }

  int total_concat_length;
  string total_concat;
  {
    vector<int> each_concat_length(comm_size);
    int local_concat_length = local_concat.size();
    if (center_node < 0)
      MPI_Allgather(&local_concat_length, 1, MPI_INT, each_concat_length.data(), 1, MPI_INT, comm);
    else
      MPI_Gather(&local_concat_length, 1, MPI_INT, each_concat_length.data(), 1, MPI_INT, center_node, comm);
    vector<int> displs(comm_size);
    displs[0] = 0;
    for (int i = 1; i < displs.size(); ++i) 
      displs[i] = displs[i - 1] + each_concat_length[i - 1];

    total_concat_length = displs.back() + each_concat_length.back();
    total_concat.resize(total_concat_length);
    if (center_node < 0)
      MPI_Allgatherv(
          (void*)local_concat.data(), local_concat_length, MPI_CHAR,
	  (void*)total_concat.data(), each_concat_length.data(), displs.data(), MPI_CHAR,
	  comm
        );
    else
      MPI_Gatherv(
          (void*)local_concat.data(), local_concat_length, MPI_CHAR,
	  (void*)total_concat.data(), each_concat_length.data(), displs.data(), MPI_CHAR,
	  center_node, comm
        );
  }

  if ((center_node < 0) || (center_node == comm_rank)) {
    all_entities.clear();
    int pos = 0;
    for (int i = 0; i < total_num; ++i) {
      string buf = total_concat.substr(pos, all_length[i]);
      all_entities.push_back(buf);
      pos += all_length[i];
    }
  }
}

void all2all_protobuf(
    std::vector<std::string> &send_entities, 
    std::vector<std::string> &recv_entities,
    MPI_Comm comm
  )
{
  int comm_size;
  int comm_rank;
  MPI_Comm_size(comm, &comm_size);
  MPI_Comm_rank(comm, &comm_rank);

  // send and receive protobuf
  std::vector<int> send_length;
  std::string send_concat;
  for (int i = 0; i < comm_size; i ++) {
    send_length.push_back(send_entities[i].size());
    send_concat += send_entities[i];
  }

  std::vector<int> recv_length(comm_size);
  MPI_Alltoall((void*)send_length.data(), 1, MPI_INT, (void*)recv_length.data(), 1, MPI_INT, comm);

  std::vector<int> send_displs(comm_size), recv_displs(comm_size);
  send_displs[0] = recv_displs[0] = 0;
  for (int i = 1; i < comm_size; i ++) {
    send_displs[i] = send_displs[i-1] + send_length[i-1];
    recv_displs[i] = recv_displs[i-1] + recv_length[i-1];
  }

  std::string recv_concat;
  recv_concat.resize(recv_displs.back() + recv_length.back());
  MPI_Alltoallv((void*)send_concat.data(), send_length.data(), send_displs.data(), MPI_CHAR, 
    (void*)recv_concat.data(), recv_length.data(), recv_displs.data(), MPI_CHAR, comm);

  // collect results
  recv_entities.clear();
  int pos = 0;
  for (int i = 0; i < comm_size; i ++) {
    if (recv_length[i] == 0) continue;
    std::string buf = recv_concat.substr(pos, recv_length[i]);
    recv_entities.push_back(buf);

    pos += recv_length[i];
  }
}  
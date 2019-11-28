#include <cmath>
#include <utility>

#include "histogram.h"

struct FileHeader {
  int nranks; // number of ranks
  int nstats; // number of stats per rank
  int nfloats; // number of 64bit floating point numbers per stat
};

struct Stats {
  double nsamples; // use double for simplicity
  double mean; // arithmetic mean
  double sdev; // standard deviation
  double accu; // accumulated value
  double vmin; // minimum value
  double vmax; // maximum value
};

static std::map<std::string, Stats> ComputeStats(
    const std::map<std::string, std::vector<double>> &samples) {
  std::map<std::string, Stats> stats;
  for (auto &key : samples) {
    Stats stat = {};
    const std::vector<double> &data = key.second;
    stat.nsamples = data.size();
    if (data.size() < 2) {
      stats[key.first] = stat;
      continue;
    }
    double sum = 0.0;
    double vmin = data[0];
    double vmax = data[0];
    for (auto s : data) {
      sum += s;
      vmin = (s < vmin) ? s : vmin;
      vmax = (s > vmax) ? s : vmax;
    }
    const double mean = sum / data.size();

    double sdev = 0.0;
    for (auto s : data) {
      const double zmean = s - mean;
      sdev += zmean * zmean;
    }
    sdev = std::sqrt(sdev / (data.size() - 1));

    stat.mean = mean;
    stat.sdev = sdev;
    stat.accu = sum;
    stat.vmin = vmin;
    stat.vmax = vmax;
    stats[key.first] = stat;
  }
  return stats;
}

void Histogram::Consolidate() {
  int rank, size;
  MPI_Comm_rank(comm_, &rank);
  MPI_Comm_size(comm_, &size);

  std::vector<size_t> char_bytes;
  std::vector<std::string> keys;
  size_t sum_char = 0;
  for (auto &key : samples_) {
    const size_t s = key.first.length();
    sum_char += s;
    char_bytes.push_back(s);
    keys.push_back(key.first);
  }

  const size_t nvalues = sizeof(Stats) / sizeof(double);
  const auto stats = ComputeStats(samples_);

  FileHeader header;
  header.nranks = size;
  header.nstats = static_cast<int>(stats.size());
  header.nfloats = nvalues;

  const size_t data_start =
      sizeof(FileHeader) + char_bytes.size() * sizeof(int) + sum_char;
  const size_t my_offset =
      data_start + rank * stats.size() * nvalues * sizeof(double);

  MPI_File fh;
  MPI_Status st;
  std::string name = "hist_" + name_ + ".bin";
  if (MPI_File_open(
          comm_, name.c_str(),
          MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN,
          MPI_INFO_NULL, &fh)) {
    throw std::runtime_error("Consolidate: Can not open MPI file");
  }

  if (0 == rank) {
    // header and keys
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    MPI_File_write(fh, &header, sizeof(FileHeader), MPI_BYTE, &st);
    for (size_t i = 0; i < samples_.size(); ++i) {
      int cb = static_cast<int>(char_bytes[i]);
      MPI_File_write(fh, &cb, 1, MPI_INT, &st);
      MPI_File_write(fh, keys[i].c_str(), cb, MPI_CHAR, &st);
    }
  } else {
    MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
  }
  for (const auto &s : stats) {
    // data
    MPI_File_write(fh, &s.second, sizeof(Stats), MPI_BYTE, &st);
  }
  MPI_File_close(&fh);
}

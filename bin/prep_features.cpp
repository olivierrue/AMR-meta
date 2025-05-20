#include <iostream>
#include <string>
#include <unordered_map>
#include <fstream> 
#include <vector>
#include <algorithm>

int main (int argc, char** argv)
{
  std::vector<std::string> ref_kmers;
  std::vector<std::string> ref_kmers_rc;
  std::string single_kmer;
  std::ifstream file_kmers(argv[4]);
  std::ifstream file_kmers_rc(argv[5]);

  while (std::getline(file_kmers, single_kmer)) {
    ref_kmers.push_back(single_kmer);
  }

  while (std::getline(file_kmers_rc, single_kmer)) {
    ref_kmers_rc.push_back(single_kmer);
  }

  std::unordered_map<std::string, uint> kmer_map;
  for (int i = 0; i < ref_kmers.size(); i++) {
    kmer_map[ref_kmers.at(i)] = i;
    kmer_map[ref_kmers_rc.at(i)] = i;
  }

  std::ifstream file_1(argv[1]);
  std::ifstream file_2(argv[2]);
  std::string str_1, str_2;
  std::string full_str;
  const unsigned int k = 13;
  const unsigned int l_ref = 138260;
  unsigned int l = 2;
  const unsigned int max_pairs_per_file = 500000;
  unsigned int file_index = 0;

  std::string readname;

  while (std::getline(file_1, str_1)) {
    std::getline(file_2, str_2);
    l++;

    if ((l - 2) % 4 == 0) { // ligne 1 d’un bloc FASTQ (nom du read)
      readname = str_1[0] == '@' ? str_1.substr(1) : str_1;
    }

    if (l % 4 == 0) {
      std::transform(str_1.begin(), str_1.end(), str_1.begin(), ::toupper); 
      std::transform(str_2.begin(), str_2.end(), str_2.begin(), ::toupper); 

      bool kmer_arr[l_ref] = {false};
      for (int i = 0; i < (str_1.size() - k + 1); i++) {
        auto sub = str_1.substr(i, k);
        if (kmer_map.count(sub) > 0) {
          kmer_arr[kmer_map[sub]] = true;
        }
      }

      for (int i = 0; i < (str_2.size() - k + 1); i++) {
        auto sub = str_2.substr(i, k);
        if (kmer_map.count(sub) > 0) {
          kmer_arr[kmer_map[sub]] = true;
        }
      }

      std::string kmer_str;
      bool celo = false;
      for (int i = 0; i < l_ref; i++) {
        if (kmer_arr[i]) {
          kmer_str += std::to_string(i + 1) + ",";
          celo = true;
        }
      }

      if (celo) {
        kmer_str.insert(0, readname + ","); // identifiant du read en première colonne
        kmer_str += '\n';
        full_str += kmer_str;
      }
    }

    if (l % (4 * max_pairs_per_file) == 0) {
      std::ofstream outfile(std::string(argv[3]) + "_" + std::to_string(file_index) + ".csv");
      if (outfile.is_open()) {
        full_str += "END\n"; // peut servir pour marquer la fin si besoin
        outfile << full_str;
        outfile.close();
        file_index++;
        full_str.clear();
      }
    }
  }

  std::ofstream outfile(std::string(argv[3]) + "_" + std::to_string(file_index) + ".csv");
  if (outfile.is_open()) {
    full_str += "END\n";
    outfile << full_str;
    outfile.close();
  }

  return 0;
}

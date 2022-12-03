#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include "tableio.h"

const std::string dirname = "fqtables/";
const std::string TEST_TABLE = "test_table";


/////////////////
// Saving

void write_table(unsigned* table, int size, std::string fname){

  // Write data to file
  std::ofstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size; i++)
    file.write((char*)(&table[i]), sizeof(unsigned));
  
  file.close();
  return;
}

void write_table(int* table, int size, std::string fname){

  // Write data to file
  std::ofstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size; i++)
    file.write((char*)(&table[i]), sizeof(int));
  
  file.close();
  return;
}


void write_table(unsigned** table, int size1, int size2, std::string fname){

  // Write data to file
  std::ofstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size1; i++)
    for (unsigned j = 0; j < size2; j++)
      file.write((char*)(&table[i][j]), sizeof(unsigned));
  
  file.close();
  return;
}

void write_table(unsigned*** table, int size1, int size2, int size3, std::string fname){

  // Write data to file
  std::ofstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size1; i++)
    for (unsigned j = 0; j < size2; j++)
      for (unsigned k = 0; k < size3; k++)
	file.write((char*)(&table[i][j][k]), sizeof(unsigned));
  
  file.close();
  return;
}

/////////////////
// Loading
//

unsigned* read_table(int size, std::string fname){

  // malloc.
  unsigned* table = new unsigned[size];
  for (unsigned i = 0; i < size; i++)
    table[i] = 0;

  // Extract data from file.
  std::ifstream file;
  file.open(dirname + fname, std::ios_base::binary);
  file.read((char*)(table), sizeof(unsigned) * size);
  //for (unsigned i = 0; i < size; i++)
  //  file.read((char*)(&table[i]), sizeof(unsigned));
  
  file.close();
  return table;
}

int* read_table(int size, std::string fname, int why){

  // malloc.
  int* table = new int[size];
  for (int i = 0; i < size; i++)
    table[i] = 0;
  
  // Extract data from file
  std::ifstream file;
  file.open(dirname + fname, std::ios_base::binary);
  file.read((char*)(table), sizeof(int) * size);
  //for (unsigned i = 0; i < size; i++)
  
  file.close();
  return table;
}


unsigned** read_table(int size1, int size2, std::string fname){

  unsigned** table = new unsigned*[size1];
  for (unsigned i = 0; i < size1; i++) {
    table[i] = new unsigned[size2];
    for (unsigned j = 0; j < size2; j++)
      table[i][j] = 0;
  }
  
  // Extract data from file
  std::ifstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size1; i++)
    file.read((char*)(table[i]), sizeof(unsigned) * size2);
    // for (unsigned j = 0; j < size2; j++)
    //   file.read((char*)(&table[i][j]), sizeof(unsigned));
  
  file.close();
  return table;
}


unsigned*** read_table(int size1, int size2, int size3, std::string fname){

  // Extract data from file.
  std::ifstream file;
  int N = size1 * size2 * size3;
  unsigned* blah = read_table(N, fname);
  file.close();
  
  // Set up the pointer table.
  unsigned*** table = new unsigned**[size1];
  int M = 0;
  for (unsigned i = 0; i < size1; i++) {
    table[i] = new unsigned*[size2];
    for (unsigned j = 0; j < size2; j++) {
      table[i][j] = (unsigned *)(blah + M); // This changes when reversing indexing.
      M += size3;
    }
  }
  return table;
}


/////////////////
// Testing
//

void remove_test_table() {
  std::string namen = dirname + TEST_TABLE;
  std::remove(namen.c_str());
  return;
}

int basic_test_table() {

  int q = 4;
  unsigned** table = new unsigned*[q];

  for (int i = 0; i < q; i++){
    table[i] = new unsigned[q];
    for (int j = 0; j < q; j++){
      table[i][j] = i + j;
    }
  }

  write_table(table, q, q, TEST_TABLE);
  unsigned** new_table = read_table(q, q, TEST_TABLE);

  for (int i = 0; i < q; i++){
    for (int j = 0; j < q; j++){
      if (table[i][j] != new_table[i][j])
	std::cerr << i << " " << j << " " << i + j << " " << new_table[i][j];      
    }
  }
  
  return 0;
}


int test_1_table(unsigned* table, int size1) {
  write_table(table, size1, TEST_TABLE);
  unsigned* new_table = read_table(size1, TEST_TABLE);

  for (int i = 0; i < size1; i++)
    if (table[i] != new_table[i])
      std::cerr << i << " " << table[i] << " " << new_table[i] << std::endl;

  return 0;
}

int test_2_table(unsigned** table, int size1, int size2) {
  write_table(table, size1, size2, TEST_TABLE);
  unsigned** new_table = read_table(size1, size2, TEST_TABLE);

  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size2; j++)
      if (table[i][j] != new_table[i][j])
	std::cerr << i << " " << j << " " << table[i][j] << " " << new_table[i][j] << std::endl;

  return 0;
}

int test_3_table(unsigned*** table, int size1, int size2, int size3) {

  write_table(table, size1, size2, size3, TEST_TABLE);
  unsigned*** new_table = read_table(size1, size2, size3, TEST_TABLE);

  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size2; j++)
      for (int k = 0; k < size3; k++)
	if (table[i][j][k] != new_table[i][j][k])
	  std::cerr << i << " " << j << " " << k << " " << table[i][j][k] << " " << new_table[i][j][k] << std::endl;

  return 0;
}

int compare_3_table(unsigned*** table1, unsigned*** table2, int size1, int size2, int size3) {

  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size2; j++)
      for (int k = 0; k < size3; k++)
	if (table1[i][j][k] != table2[i][j][k])
	  std::cerr << i << " " << j << " " << k << " " << table1[i][j][k] << " " << table2[i][j][k] << std::endl;

  return 0;
}

int compare_2_table(unsigned** table1, unsigned** table2, int size1, int size2) {

  for (int i = 0; i < size1; i++)
    for (int j = 0; j < size2; j++)
      if (table1[i][j] != table2[i][j])
        std::cerr << i << " " << j << " " << table1[i][j] << " " << table2[i][j] << std::endl;

  return 0;
}

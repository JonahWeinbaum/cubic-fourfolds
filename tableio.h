void write_table(unsigned* table, int size, std::string fname);
void write_table(int* table, int size, std::string fname);
void write_table(unsigned** table, int size1, int size2, std::string fname);
void write_table(unsigned*** table, int size1, int size2, int size3, std::string fname);

int* read_table(int size, std::string fname, int why);
unsigned* read_table(int size, std::string fname);
unsigned** read_table(int size1, int size2, std::string fname);
unsigned*** read_table(int size1, int size2, int size3, std::string fname);

int basic_test_table();
int test_1_table(unsigned* table, int size1);
int test_2_table(unsigned** table, int size1, int size2);
int test_3_table(unsigned*** table, int size1, int size2, int size3);

int compare_3_table(unsigned*** table1, unsigned*** table2, int size1, int size2, int size3);

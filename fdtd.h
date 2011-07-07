void init();
double get_delta_t();
void open_output_files();
void calc_fdtd(int step, int last_step);
void close_output_files();
double eabs(int i, int j, int k);


/* header file for lib_pio.c */

void read_pdb_silent(char *, unsigned *, PDB **);
void read_pdb(char *, unsigned *, PDB **);
int coord_precision(double);
void write_pdb(char *, unsigned, PDB *);
void append_pdb(char *, unsigned, PDB *);
void fld2s(char *, char * );
void fld2i(char *, unsigned, int *);
void get_fld(char *, unsigned, unsigned, char * );
void print_space(FILE *, unsigned);



#define MAXLEN 10000
#include <string>
using namespace std;

namespace TMalignC{
class TMalign_wrapper
{
public: 

	char version[20];							   //for version name
	double D0_MIN;                             //for d0
	double Lnorm;                              //normalization length
	double score_d8, d0, d0_search, dcu0;      //for TMscore search
	double **score;            			       //Input score table for dynamic programming
	bool   **path;                             //for dynamic programming  
	double **val;                              //for dynamic programming  
	int    xlen, ylen, minlen;                 //length of proteins
	double **xa, **ya;                         //for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2]
                                           //in general, ya is regarded as native structure --> superpose xa onto ya
	int    *xresno, *yresno;                   //residue numbers, used in fragment gapless threading 
	double **xtm, **ytm;                       //for TMscore search engine
	double **xt;                               //for saving the superposed version of r_1 or xtm
	char   *seqx, *seqy;                       //for the protein sequence 
	int    *secx, *secy;                       //for the secondary structure 
	double **r1, **r2;                         //for Kabsch rotation 
	double t[3], u[3][3];                      //Kabsch translation vector and rotation matrix


	//output values. stores TMalign results that are
	//obtained by cython
	double rmsd_out, TM1_out, TM2_out, number_aligned;
	double t0_out[3], u0_out[3][3], seqid_out;
	std::string seqX_align_out, align_out, seqY_align_out; 

	//argument variables
	char out_reg[MAXLEN];
	double Lnorm_ass, Lnorm_d0, d0_scale, d0A, d0B, d0u, d0a;
	bool o_opt, a_opt, u_opt, d_opt, v_opt;
	double TM3, TM4, TM5;

	TMalign_wrapper();

	void char_test();
	void print_help(char *arg);
	void parameter_set4search(int xlen, int ylen);
	void parameter_set4final(double len);
	void parameter_set4scale(int len, double d_s);
	void load_PDB_allocate_memory(char *xname, char *yname);
	void free_memory();
	int score_fun8( double **xa, 
	                double **ya, 
	                int n_ali,
	                double d,
	                int i_ali[], 
	                double *score1,
	                int score_sum_method
	              );

	double TMscore8_search( double **xtm, 
	                        double **ytm,
	                        int Lali, 
	                        double t0[3],
	                        double u0[3][3],
	                        int simplify_step,
	                        int score_sum_method,
	                        double *Rcomm
	                       );
	// TMscore search engine
	// input:   two aligned vector sets: x, y
	//          scale parameter d0
	//          simplify_step: 1 or 40 or other integers
	//          score_sum_method: 0 for score over all pairs
	//                            8 for socre over the pairs with dist<score_d8
	//                                  
	//          
	// output:  the best rotaion matrix t0, u0 that results in highest TMscore


	double detailed_search( double **x,
	                        double **y, 
	                        int x_len, 
	                        int y_len, 
	                        int invmap0[],
	                        double t[3],
	                        double u[3][3],
	                        int simplify_step,
	                        int score_sum_method                        
                       );
	//Comprehensive TMscore search engine
	// input:   two vector sets: x, y
	//          an alignment invmap0[] between x and y
	//          simplify_step: 1 or 40 or other integers
	//          score_sum_method: 0 for score over all pairs
	//                            8 for socre over the pairs with dist<score_d8          
	// output:  the best rotaion matrix t, u that results in highest TMscore

	double get_score_fast(double **x, double **y, int x_len, int y_len, int invmap[]);
	//compute the score quickly in three iterations

	double get_initial( double **x, 
                    double **y, 
                    int x_len,
                    int y_len, 
                    int *y2x
                   );
	//perform gapless threading to find the best initial alignment
	//input: x, y, x_len, y_len
	//output: y2x0 stores the best alignment: e.g., 
	//y2x0[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0 
	//the jth element in y is aligned to a gap in x if i==-1

	void smooth(int *sec, int len);
	int sec_str(double dis13, double dis14, double dis15, double dis24, double dis25, double dis35);

	void make_sec(double **x, int len, int *sec);
	//1->coil, 2->helix, 3->turn, 4->strand

	void get_initial_ss(  double **x, 
					  double **y, 
					  int x_len,
					  int y_len, 
					  int *y2x
                      );
	//get initial alignment from secondary structure alignment
	//input: x, y, x_len, y_len
	//output: y2x stores the best alignment: e.g., 
	//y2x[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0 
	//the jth element in y is aligned to a gap in x if i==-1


	bool get_initial_local(  double **x, 
						 double **y, 
						 int x_len,
						 int y_len, 
						 int *y2x
						 );
	// get_initial5 in TMalign
	//get initial alignment of local structure superposition
	//input: x, y, x_len, y_len
	//output: y2x stores the best alignment: e.g., 
	//y2x[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0 
	//the jth element in y is aligned to a gap in x if i==-1

	void score_matrix_rmsd(  double **x, 
						 double **y, 
						 int x_len,
						 int y_len,
						 int *y2x
						 );
	//with invmap(i) calculate score(i,j) using RMSD rotation

	void score_matrix_rmsd_sec(  double **x, 
					 	double **y, 
						int x_len,
						int y_len,
						int *y2x
						);

	void get_initial_ssplus( double **x, 
						 double **y, 
						 int x_len,
						 int y_len,
						 int *y2x0,
						 int *y2x						
						 );
	//get initial alignment from secondary structure and previous alignments
	//input: x, y, x_len, y_len
	//output: y2x stores the best alignment: e.g., 
	//y2x[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0 
	//the jth element in y is aligned to a gap in x if i==-1

	void find_max_frag(double **x, int *resno, int len, int *start_max, int *end_max);

	double get_initial_fgt( double **x, 
						double **y, 
						int x_len,
						int y_len, 
						int *xresno,
						int *yresno,
						int *y2x
						);
	//perform fragment gapless threading to find the best initial alignment
	//input: x, y, x_len, y_len
	//output: y2x0 stores the best alignment: e.g., 
	//y2x0[j]=i means:
	//the jth element in y is aligned to the ith element in x if i>=0 
	//the jth element in y is aligned to a gap in x if i==-1


	double DP_iter( double **x,
	                double **y, 
	                int x_len, 
	                int y_len, 
	                double t[3],
	                double u[3][3],
			int invmap0[],
			int g1,
			int g2,
			int iteration_max                                   
			);
	//heuristic run of dynamic programing iteratively to find the best alignment
	//input: initial rotation matrix t, u
	//       vectors x and y, d0
	//output: best alignment that maximizes the TMscore, will be stored in invmap

	void output_superpose(char *xname,
					  char *yname,
					  int x_len,
					  int y_len,
					  double t[3],
					  double u[3][3],
					  double rmsd,
					  double d0_out,
					  int m1[], 
					  int m2[],
					  int n_ali8,
					  double seq_id,
					  double TM_0,
					  double Lnorm_0,
					  double d0_0
					 );

	void output_results(char *xname,
					 char *yname,
					 int x_len,
					 int y_len,
					 double t[3],
					 double u[3][3],
					 double TM1,
					 double TM2,
					 double rmsd,
					 double d0_out,
					 int m1[], 
					 int m2[],
					 int n_ali8,
					 int n_ali,
					 double TM_0,
					 double Lnorm_0,
					 double d0_0
					 );
	//output the final results

	void PrintErrorAndQuit(std::string sErrorString);

	template <class A> void NewArray(A *** array, int Narray1, int Narray2);

	template <class A> void DeleteArray(A *** array, int Narray);

	char AAmap(std::string AA);
	void AAmap3(char A, char AA[3]);
	void get_xyz(std::string line, double *x, double *y, double *z, char *resname, int *no);
	int get_PDB_len(char *filename);
	int read_PDB(char *filename, double **a, char *seq, int *resno);
	int get_ligand_len(char *filename);
	int read_ligand(char *filename, double **a, char *seq, int *resno);
	double dist(double x[3], double y[3]);
	double dot(double *a, double *b);
	void transform(double t[3], double u[3][3], double *x, double *x1);
	void do_rotation(double **x, double **x1, int len, double t[3], double u[3][3]);
	void output_align1(int *invmap0, int len);
	int output_align(int *invmap0, int len);
	void NWDP_TM(int len1, int len2, double gap_open, int j2i[]);
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
	//Input: score[1:len1, 1:len2], and gap_open
	//Output: j2i[1:len2] \in {1:len1} U {-1}
	//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical
	void NWDP_TM(double **x, double **y, int len1, int len2, double t[3], double u[3][3], double d02, double gap_open, int j2i[]);
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
	//Input: vectors x, y, rotation matrix t, u, scale factor d02, and gap_open
	//Output: j2i[1:len2] \in {1:len1} U {-1}
	//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical
	void NWDP_TM(int *secx, int *secy, int len1, int len2, double gap_open, int j2i[]);
	//NW dynamic programming for alignment
	//not a standard implementation of NW algorithm
    //Input: secondary structure secx, secy, and gap_open
    //Output: j2i[1:len2] \in {1:len1} U {-1}
	//path[0:len1, 0:len2]=1,2,3, from diagonal, horizontal, vertical
	bool Kabsch(double **x, 
            double **y, 
            int n, 
            int mode,
            double *rms,
            double t[3],
            double u[3][3]
            );
	void allocate_memory();
	void run_tmalign();
	void run_main(int argc, char *argv[]);
        void obtain_results(char *xname,
                                         char *yname,
                                         int x_len,
                                         int y_len,
                                         double t[3],
                                         double u[3][3],
                                         double TM1,
                                         double TM2,
                                         double rmsd,
                                         double d0_out,
                                         int m1[],
                                         int m2[],
                                         int n_ali8,
                                         int n_ali,
                                         double TM_0,
                                         double Lnorm_0,
                                         double d0_0
                                         );

};
}

#ifndef SBLAST_H_
#define SBLAST_H_

#include<stdio.h>
#include<string.h>
#include"svector.h"

const int MEM_BLOCK=16;
const int LOCUS_LEN=128;

const unsigned SB_INVERSION = 0x0001;
const unsigned SB_BOTH_GAP = 0x0002;
const unsigned SB_MULTI = 0x0004;
const unsigned SB_FORWARD = 0x0008; // not for end users
const unsigned SB_NOT_MASK_REPEAT = 0x0010;
const unsigned SB_BRIEF = 0x0020;
const unsigned SB_LINEAR_CLUSTER = 0x0040; // colinear clustering
const unsigned SB_CROSS_SBJCT = 0x0080;
const unsigned SB_CLUSTER = 0x0100;

const int SB_UNALIGNED = -1;
const int SB_ALIGNED = -2;

typedef int node_t;
typedef int score_t;
typedef int length_t;

struct aBasicBlock
{
	length_t qt, qp, st, sp;
	score_t score;
	bool dir;
//	double e; // e-value, but it is still useless in present version.
};
struct aBlock : public aBasicBlock
{
	int repeat;
	int s_locus_id;
};
// there is another version of "<" in flag_all.cc
inline bool operator< (const aBlock &a, const aBlock &b)
{
	if (a.s_locus_id < b.s_locus_id) return true;
	if (a.s_locus_id > b.s_locus_id) return false;
	return a.qt < b.qt;
}

struct aLine
{
	length_t num[9];
	char locus[2][LOCUS_LEN];
	double e;
};

class Blocks
{
	FILE *fp;
	node_t max;
	length_t max_send;
	aLine *last,*curr;
	void read_line(aLine*);
	void add_first();
	void add(int*, double);
protected:
	node_t num; // number of blocks
	char *q_locus;
	length_t q_len;
	svector<char*> s_locus;
	svector<length_t> s_len;
	aBlock *blocks;

	virtual void reset();
public:
	Blocks(FILE *f=stdin);
	virtual ~Blocks(void);
	node_t read_blocks(void);
	inline int get_num(void) { return(num); };
};

struct VexEdge
{
	node_t num, max, *array;
	length_t score;
	inline void add(node_t n);
	inline void rewind(void) { num=0; };
};

struct SBResultStruct : public aBasicBlock
{
	node_t num, max;
	aBlock **list;
};

class SortBlast:public Blocks
{
protected:
	node_t ve_max; // maximum size of ve;
	svector<same_pair<length_t> > rep_regions; // repeated regions, generated in flag_repeat()
	svector<aBlock*> rbs;  // Repeated BlockS
	svector<aBlock*> wbs;  // Working BlockS
	
	node_t w_num; // number of working blocks
	VexEdge *ve; // vertics and related edges
	
	// in flag_rep.cc
	void flag_repeat(); // flag repeats
	void reset(); // virtual function

	// This function has been revised to give a older output format. This format will be larger and do
	// not contain repeat information, but it is more friendly and comfortable to eyes.
	// old=-1 means new format.
	void output(const SBResultStruct &, FILE *fpout, int); // output the line A ......
	bool write_result(node_t*, node_t, SBResultStruct*);
	
	// in sortBlast.cc
	inline bool isconsist(aBlock*, aBlock*); // judge if the two blocks are conflict with each other
	int optimal(node_t *); // dynamic programming
	void build_list(void); // build graph for optimization
	int back_trace(node_t *,node_t); // back trace the optimal path
	double penalty(aBlock*, aBlock*); // penalty function
	score_t cal_score(node_t *stack, node_t size); // get rid of the repeated scoring at overlapped region
	node_t get_optimal(node_t *stack); // interface for optimal() build_list() penalty() and back_trace()
public:
	SortBlast(FILE *f, unsigned t);
	~SortBlast(void);
};

extern double SB_overlap_ratio;
extern double SB_long_penalty;
extern int SB_overlap_len;
extern int SB_score_min_len;
extern int SB_max_gap_len;
extern bool SB_half_gap;
extern int SB_rep_bound; // in flag_rep.cc
extern int SB_rep_depth; // in flag_rep.cc
extern int SB_max_cluster_gap; // in cluster.cc
extern int SB_min_chain_score; // write_result.cc

extern unsigned sb_flag; // in readFile.cc

inline bool operator<(const aBlock &a,const aBlock &b);
score_t cal_part_score(aBlock*, aBlock*);

#endif

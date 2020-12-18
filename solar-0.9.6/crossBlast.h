#ifndef CROSSBLAST_H_
#define CROSSBLAST_H_

#include <stdio.h>
#include"sBlast.h"
#include"svector.h"

class CrossBlast:public SortBlast
{
protected:
	inline bool isconflict(aBlock*, aBlock*);
	bool gen_work_blocks(node_t*, node_t);
	int gen_first_work_blocks(int*);
	void gen_alignment(FILE*, SBResultStruct*, int);
public:
	void handle(FILE*);
	CrossBlast(FILE *fp, unsigned);
	~CrossBlast(void);
};

void sb_cluster();

#endif

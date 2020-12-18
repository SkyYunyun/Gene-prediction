#include<unistd.h>
#include<stdio.h>
#include<stdlib.h>
#include"crossBlast.h"

void handle(FILE *fp, unsigned t, FILE *fpout)
{
	CrossBlast sb(fp, t);
	sb.handle(fpout);
}
void usage(void)
{
	fprintf(stderr,"\n");
	fprintf(stderr,"Program : solar (Sorting Out Local Alignment Result)\n");
	fprintf(stderr,"Version : 0.9.6, on 08 November, 2006\n");
	fprintf(stderr,"Contact : liheng@genomics.org.cn\n\n");
	fprintf(stderr,"Usage   : snap [options] [input [output]]\n\n");
	fprintf(stderr,"Options : -n INUM   maximum gap length, default is %d\n", SB_max_gap_len);
	fprintf(stderr,"          -c        cluster and construct multi-blocks\n");
	fprintf(stderr,"          -C        not examine the overlap in query (may be VERY SLOW)\n");
	fprintf(stderr,"          -h        help\n\n");
	fprintf(stderr,"Advanced: -s INUM   maximum cluster gap, default is %d\n", SB_max_cluster_gap);
	fprintf(stderr,"          -m INUM   minimum score of a chain, default is %d\n", SB_min_chain_score);
	fprintf(stderr,"          -b        detailed output format\n");
	fprintf(stderr,"          -p FNUM   penalty for long gap, default is %.2f\n",SB_long_penalty);
	fprintf(stderr,"          -g INUM   minimum gap length for penalty, default is %d\n",SB_score_min_len);
	fprintf(stderr,"                    penalty = long_gap_pen * log(gap_len / min_gap_len)\n");
	fprintf(stderr,"          -u INUM   maximum size flanking the repeat boundary, default is %d\n",
			SB_rep_bound);
	fprintf(stderr,"          -d INUM   minimum depth for repeats (-1 stands for no masking),\n");
	fprintf(stderr,"                    default is %d\n", SB_rep_depth);
	fprintf(stderr,"          -r FNUM   threshold for overlap ratio, default is %.2f\n",SB_overlap_ratio);
	fprintf(stderr,"          -l INUM   maximum length of overlap length, -1 for no limit\n");
	fprintf(stderr,"                    default is -1\n");
	fprintf(stderr,"          -v        permitting inversion in one block\n");
	fprintf(stderr,"          -t        use the difference of the gaps\n\n");
	fprintf(stderr,"Comment : Input file format:\n");
	fprintf(stderr,"            Qname Qlen Sname Slen Qstart Qstop Sstart Sstop score e-value\n\n");
	fprintf(stderr,"            At present, e-value is useless. Qlen and Slen are not necessary,\n");
	fprintf(stderr,"          either. These three fields can be safely set as zero. Then Qlen\n");
	fprintf(stderr,"          and Slen will be estimated by snap in output, not precisely though.\n\n");
	fprintf(stderr,"          Concise output format:\n");
	fprintf(stderr,"            Qname Qlen Qstart Qstop strand Sname Slen Sstart Sstop #blocks \\\n");
	fprintf(stderr,"               total_score Qstart,Qstop;...; Sstart,Sstop;...; score;...;\n\n");
	fprintf(stderr,"          Detailed output format:\n");
	fprintf(stderr,"            Q	Qname Qlen\n");
	fprintf(stderr,"            R	#repeats Rstart,Rend;...;\n");
	fprintf(stderr,"            S	Sname Slen\n");
	fprintf(stderr,"            C#	Qstart Qend * Sstart Send #blocks total_score\n");
	fprintf(stderr,"            A	Qstart Qend strand Sstart Send #blocks total_score \\\n");
	fprintf(stderr,"                Qstart,Qend;...; Sstart,Send;...; score;...;\n");
	fprintf(stderr,"            T	#alignments total_score\n");
	fprintf(stderr,"            F	#remained Qstart,Qend;...; Tstart,Tend;...; score;...;\n\n");
	fprintf(stderr,"          where Q stands for Query, R for Repeat, S Subject, C Cluster,\n");
	fprintf(stderr,"          A Alignment, T Total and F for Fragment.\n\n");
	exit(1);
}

int main(int argc,char *argv[])
{
	int c;
	FILE *fp = stdin, *fpout = stdout;
	sb_flag = SB_BRIEF; // default flag
	while ((c = getopt(argc,argv,"n:s:m:cCbhp:g:u:d:r:l:vt")) >= 0) {
		switch(c) {
			case 'n': SB_max_gap_len = atoi(optarg); break;
			case 's': SB_max_cluster_gap = atoi(optarg); break;
			case 'm': SB_min_chain_score = atoi(optarg); break;
			case 'c': sb_flag |= SB_MULTI | SB_CLUSTER; break;
			case 'C': sb_flag |= SB_MULTI | SB_CROSS_SBJCT; sb_flag |= SB_CLUSTER; break;
			case 'b': sb_flag &= ~SB_BRIEF; break;
			case 'p': SB_long_penalty = atof(optarg); break;
			case 'g': SB_score_min_len = atoi(optarg); break;
			case 'u': SB_rep_bound = atoi(optarg); break;
			case 'd': SB_rep_depth = atoi(optarg);
					  if (SB_rep_depth < 0) sb_flag |= SB_NOT_MASK_REPEAT;
					  break;
			case 'r': SB_overlap_ratio = atof(optarg); break;
			case 'l': SB_overlap_len = atoi(optarg);
					  if (SB_overlap_len < 0) SB_overlap_len = 1<<30;
					  break;
			case 'v': sb_flag |= SB_INVERSION; break;
			case 't': sb_flag |= SB_BOTH_GAP; break;
			case 'h':usage();
		}
	}
	if(optind<argc) {
		fp=fopen(argv[optind],"r");
		if(!fp) {
			fprintf(stderr,"Cannot open the file %s\n",argv[optind]);
			exit(1);
		}
		optind++;
		if (optind < argc) {
			fpout = fopen(argv[optind], "w+");
			if (!fpout) {
				fprintf(stderr, "Cannot create the file %s\n",argv[optind]);
				exit(1);
			}
		}
	}
	handle(fp, sb_flag, fpout);
	fclose(fp);
	fclose(fpout);
	return(0);
}

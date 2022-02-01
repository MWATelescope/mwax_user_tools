/*
Read MWAX subfile(s) and reinterleave the data into a stream.
Remove the header and metadata blocks and optionally output them to different files
Randall Wayth. Dec 2021.

The output data stream has time varying most slowly
*/

/*
MWAX subfiles are comprised of:
4096 byte header (text)
161 blocks, where each block is of size blocksize=n_inputs*n_samplesperblock*n_bytespersample
where
n_inputs=256 # is actually 128 inputs * 2 pols
n_samplesperblock=64000
n_bytespersample=2 # 1 byte each for I and Q

The first block is reserved for metadata, followed by 160 blocks of I/Q samples.
The data ordering in a block is (slowest->fastest changing) [n_times][n_inputs][n_pols][n_samples]
so that there are two time indices. n_samples is 64000 samples all from the same input.
We want to re-arrange this so that the output ordering is
[time][input] (or equivalently [time][input][pol])
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#define DEFAULT_SAMP_SIZE         2
#define DEFAULT_SAMP_PERBLOCK     64000
#define DEFAULT_NINP              256
#define HEADER_SIZE               4096

/* globals */
int bytes_per_samp=DEFAULT_SAMP_SIZE;   
int n_inputs=DEFAULT_NINP;       // including 2 pols, but we don't need to explicitly worry about them
int n_samp_per_block=DEFAULT_SAMP_PERBLOCK;
int n_blocks_per_subfile=160;
int debug=0;

char *infilename=NULL, *outfilename=NULL, *outfilename_meta=NULL,*outfilename_header=NULL;
FILE *fpd,*fpin=NULL, *fpout_data=NULL, *fpout_header=NULL, *fpout_meta=NULL;


void print_usage(char * const argv[]) {
    fprintf(stderr,"Usage:\n%s [options] \n",argv[0]);
    fprintf(stderr,"\t-s num        \tspecify sample size in bytes. default: %d\n",DEFAULT_SAMP_SIZE);
    fprintf(stderr,"\t-i filename\tinput file name. default: stdin\n");
    fprintf(stderr,"\t-o filename\toutput file name. default: stdout\n");
    fprintf(stderr,"\t-h filename\toutput file for the headers. Default: don't write\n");
    fprintf(stderr,"\t-m filename\toutput file for the metadata blocks. Default: don't write\n");
    fprintf(stderr,"\t-s number\tNumber of samples per block. Default: %d\n",DEFAULT_SAMP_PERBLOCK);
    fprintf(stderr,"\t-n number\tTotal number of inputs (antennas*pols). Default: %d\n",DEFAULT_NINP);
    fprintf(stderr,"\t-d            \twrite debug and runtime info to stderr\n");
    exit(0);
}

void parse_options(int argc, char * const argv[]) {
    const char optstring[]="do:s:i:h:m:n:";
    int c;
    
    while ((c=getopt(argc,argv,optstring)) != -1) {
        switch(c) {
            case 's':
                n_inputs = atoi(optarg);
                if (n_inputs <=0 || n_inputs > 1024) {
                    fprintf(stderr,"bad num inputs: %d\n",n_inputs);
                    print_usage(argv);
                }
                break;
            case 'n':
                n_samp_per_block = atoi(optarg);
                if (n_samp_per_block <=0 || n_samp_per_block > 1000000) {
                    fprintf(stderr,"bad num samps per block: %d\n",n_samp_per_block);
                    print_usage(argv);
                }
                break;
            case 'd':
                debug += 1;
                break;
            case 'o':
                outfilename=optarg;
                break;
            case 'i':
                infilename=optarg;
                break;
            case 'm':
                outfilename_meta=optarg;
                break;
            case 'h':
                outfilename_header=optarg;
                break;
            default:
                fprintf(stderr,"unknown option %c\n",c);
                print_usage(argv);
        }
    }
    /* process remaining options */
   if (debug) {
        fprintf(stderr,"N inputs: %d\n",n_inputs);
        fprintf(stderr,"Samples per block: %d\n",n_samp_per_block);
    }
}

/* 
   Reorder/transpose samples within a single sub-block of n_samp_per_block
   we want to copy chunks of samples from the input to the output, doing a corner turn in the process
   Do strided reads and coalesced writes, assuming this is the most efficient approach
   input  index = samp + inp*nsamp + time*nsamp*ninp
   output index = inp + samp*ninp + time*nsamp*ninp
*/
/* this version for 8+8 real/imag samples */
void do_reorder_8p8(uint16_t *in, uint16_t *out) {
    int inp,samp,inp_ind,out_ind;
    
    for (samp=0; samp<n_samp_per_block; samp++) {
        for (inp=0; inp < n_inputs; inp++) {
            inp_ind = samp + inp*n_samp_per_block;
            out_ind = inp + samp*n_inputs;
            out[out_ind] = in[inp_ind];
        }
    }
}


int main(int argc, char * const argv[]) {
    int i,done=0,n;
    char *outdata=NULL,*datablock=NULL,*header=NULL,*metablock=NULL;
    int32_t blocksize;

    fpd=stderr;

    if (argc <2) print_usage(argv);
    parse_options(argc,argv);

    blocksize=n_inputs*n_samp_per_block*bytes_per_samp;
    if(debug) {
        fprintf(stderr,"Block size %d bytes\n",blocksize);
    }

    /* open files */
    if (outfilename==NULL) {
        fpout_data=stdout;
    }
    else {
        fpout_data=fopen(outfilename,"w");
        if(fpout_data==NULL) {
            fprintf(stderr,"Cannot open output file %s\n",outfilename);
            return -1;
        }
    }

    if (infilename==NULL) {
        fpin=stdin;
    }
    else {
        fpin=fopen(infilename,"r");
        if(fpin==NULL) {
            fprintf(stderr,"Cannot open in file %s\n",infilename);
            return -1;
        }
    }

    if (outfilename_header !=NULL) {
        fpout_header = fopen(outfilename_header,"w");
        if (fpout_header==NULL) {
            fprintf(stderr,"ERROR: cannot open output file for header info %s\n",outfilename_header);
            return -1;
        }
    }
    
    if (outfilename_meta != NULL) {
        fpout_meta = fopen(outfilename_meta,"w");
        if (fpout_meta==NULL) {
            fprintf(stderr,"ERROR: cannot open output file for metadata  %s\n",outfilename_meta);
            return -1;
        }        
    }

    /* make space for input and output chunks */
    header=malloc(HEADER_SIZE);
    assert(header!=NULL);
    metablock=malloc(blocksize);
    assert(metablock!=NULL);
    datablock=malloc(blocksize);
    assert(datablock!=NULL);
    outdata=malloc(blocksize);
    assert(outdata!=NULL);

    /* do it */
    while (!done) {
        /* read header */
        n=fread(header,HEADER_SIZE,1,fpin);
        if (n==0) {
            // no more data, normal exit
            done=1;
            continue;
        }
        if (debug) {
            fprintf(fpd,"Read header\n");
        }
        /* read metablock */
        n=fread(metablock,blocksize,1,fpin);
        if (n==0) {
            fprintf(stderr,"ERROR: EOF trying to read metablock\n");
            done=1;
            continue;
        }

        /* write the header and metadata blocks if required */
        if (fpout_header != NULL) {
            fwrite(header,HEADER_SIZE,1,fpout_header);
        }
        if (fpout_meta != NULL) {
            fwrite(metablock,blocksize,1,fpout_meta);
        }

        /* loop through each time block */
        for (i=0; i<n_blocks_per_subfile; i++) {
            
            n=fread(datablock,blocksize,1,fpin);
            if(n==0) {
                fprintf(stderr,"ERROR: EOF trying to read data block %d\n",i);
            }
            /* do the transpose */
            do_reorder_8p8((uint16_t *)datablock,(uint16_t *)outdata);
            
            /* write the outpt data */
            n=fwrite(outdata,blocksize,1,fpout_data);
            if (n==0) {
                fprintf(stderr,"ERROR trying to write output for block index %d\n",i);
                done=1;
            }
        }
    }
    
    /* clean up */
    if (header!=NULL) free(header);
    if (metablock!=NULL) free(metablock);
    if (datablock!=NULL) free(datablock);
    if (outdata!=NULL) free(outdata);
    if (fpout_header !=NULL) fclose(fpout_header);
    if (fpout_meta !=NULL) fclose(fpout_meta);
    fclose(fpout_data);
    fclose(fpin);
}



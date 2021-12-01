// 
// *Very* basic utility to replace an existing subfile header with a header based on a supplied text file.
//
// To compile: 
//   gcc -Wall -Ofast mwax_update_subfile_header.c -oupdate_mwax_subfile_header -lrt
//
// Usage:
//   mwax_update_subfile_header HEADERFILE SUBFILE
//
// Where:
//   HEADERFILE is an ASCII text file containing the PSRDADA header you want for the SUBFILE.
//
// See the following link for a description of an MWAX Subfile PSRDADA header: https://wiki.mwatelescope.org/display/MP/MWAX+PSRDADA+header
// 
// TODO: make more user friendly with standard command line arg parsing and --help for usage info
// TODO: better error handling
// TODO: provide a method (or another tool) which will just take a subfile and correlator params and just update those. This will be useful for folks using the MWAX Offline correlator especially!
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <signal.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <time.h>

#define FALSE 0
#define TRUE !(FALSE)

#define BUILD 1

#define BOOL int
#define INT8 char
#define UINT8 unsigned char
#define INT16 short
#define UINT16 unsigned short
#define INT32 int
#define UINT32 unsigned int
#define INT64 long long int
#define UINT64 unsigned long long int

#define HEADER_LEN 4096

int main(int argc, char **argv)
{
//    char filename[300];

    int file_descr;

    int input_read;
    int input_write;

    char header_buffer[HEADER_LEN]={0x00};

    file_descr = open( argv[1], O_RDONLY );
    if (file_descr < 0) {
        printf( "Error opening header file:%s\n", argv[1] );
        exit(0);
    }


    input_read = read ( file_descr, header_buffer, HEADER_LEN );
    if (input_read < 20) {
        printf( "insufficient header read from file:%s.  Read returned %d.\n", argv[1], input_read );
        exit(0);
    }

    printf( "Read a %d length header from file:%s\n", input_read, argv[1] );

    close ( file_descr );

    file_descr = open( argv[2], O_WRONLY );
    if (file_descr < 0) {
        printf( "Error opening destination file:%s\n", argv[2] );
        exit(0);
    }

    input_write = write ( file_descr, header_buffer, HEADER_LEN );
    if (input_write != HEADER_LEN) {
        printf( "Write to %s failed.  Returned a value of %d.\n", argv[2], input_write );
        exit(0);
    }

    close ( file_descr );

    printf( "Header updated successfully\n" );

    exit(0);
}

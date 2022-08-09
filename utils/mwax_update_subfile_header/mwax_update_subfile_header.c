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
#include <ascii_header.h> // Supplied by PSRDADA

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

#define MAX_KV_SIZE 256

// A basic linked list struct
#define MODIFY_KV 0 /* Either modify or create keyword */
#define DELETE_KV 1 /* Delete keyword */
struct keyval
{
    char key[MAX_KV_SIZE];
    char val[MAX_KV_SIZE];
    int instruction; // MODIFY_KV, DELETE_KV
    struct keyval *next;
};

int append_keyval(struct keyval *kv, const char *key, const char *val, int instruction)
{
    if (kv == NULL) // If kv is an empty list...
    {
        kv = (struct keyval *)malloc(sizeof(struct keyval));
    }
    else
    {
        // Move to the end of the list
        while (kv->next != NULL)
            kv = kv->next;

        // Create a new node and move to it
        kv->next = (struct keyval *)malloc(sizeof(struct keyval));
        kv = kv->next;
    }

    // Update the values
    if (key != NULL)
        snprintf(kv->key, MAX_KV_SIZE, "%s", key);

    if (val != NULL)
        snprintf(kv->val, MAX_KV_SIZE, "%s", val);

    kv->instruction = instruction;
}

void usage(FILE *f, char **argv)
{
    fprintf(f, "usage: %s [-h] [-s KEY=VAL [-s ...]] [-d KEY [-d ...]] SUBFILE [SUBFILE ...]\n", argv[0]);
    fprintf(f, "\t-d KEY        Deletes KEY from header.\n");
    fprintf(f, "\t-h            Prints this help and exits.\n");
    fprintf(f, "\t-s KEY=VAL    Sets the value of KEY to VAL. If KEY does not exist, it is created.\n");
    fprintf(f, "If no -d or -s is given, then the subfile headers are printed out to stdout\n");
}

struct keyval *parse_cmdline(int argc, char **argv, int *first_filename_idx);

int main(int argc, char **argv)
{
    // Parse the command line
    int first_filename_idx;
    struct keyval *kv = parse_cmdline(argc, argv, &first_filename_idx);

    // Loop over the remaining arguments, open the files, and view/modify their headers
    int i;
    for (i = first_filename_idx; i < argc; i++)
    {
        int file_descr;

        int input_read;
        int input_write;

        char header_buffer[HEADER_LEN];

        file_descr = open( argv[i], O_RDONLY );
        if (file_descr < 0)
        {
            printf( "ERROR: Error opening header file: %s\n", argv[i] );
            exit(EXIT_FAILURE);
        }

        input_read = read(file_descr, header_buffer, HEADER_LEN);
        if (input_read < HEADER_LEN)
        {
            fprintf(stderr, "ERROR: insufficient header read from file:%s.  Read returned %d.\n", argv[i], input_read );
            exit(EXIT_FAILURE);
        }

        // If there are no -d's or -s's, then kv will be NULL.
        // In this case, just read in the headers and print them to stdout
        if (kv == NULL)
        {
            printf("%s:\n%s\n", argv[i], header_buffer);
            close ( file_descr );
            continue;
        }

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

        exit(EXIT_SUCCESS);
    }
}

/**
 * Parse the command line
 * @param[in] argc Same as main()
 * @param[in] argv Same as main()
 * @param[out] first_filename_idx The index into argv of the first filename
 * @return A pointer to a newly constructed linked list
 */
struct keyval *parse_cmdline(int argc, char **argv, int *first_filename_idx)
{
    struct keyval *kv = NULL;
    char key[MAX_KV_SIZE];
    char val[MAX_KV_SIZE];
    int n; // Used for counting scanf'd items

    int opt;
    while ((opt = getopt(argc, argv, "d:hs:")) != -1)
    {
        switch (opt)
        {
            case 'd':
                append_keyval(kv, key, NULL, DELETE_KV);
                break;
            case 'h':
                printf("Modifies subfile headers in-place\n");
                usage(stdout, argv);
                exit(EXIT_SUCCESS);
                break;
            case 's':
                n = sscanf(optarg, "%[^=]=%s", key, val);
                if (n != 2)
                {
                    printf("n=%d\n", n);
                    fprintf(stderr, "ERROR: cannot parse \"%s\" as KEY=VAL\n", optarg);
                    usage(stderr, argv);
                    exit(EXIT_FAILURE);
                }
                append_keyval(kv, key, val, MODIFY_KV);
                break;
            default:
                fprintf(stderr, "ERROR: Unrecognised option '-%c'\n", opt);
                usage(stderr, argv);
                exit(EXIT_FAILURE);
        }
    }

    // Check that there is at least one non-option argument given
    if (optind == argc)
    {
        fprintf(stderr, "ERROR: No subfiles given\n");
        usage(stderr, argv);
        exit(EXIT_FAILURE);
    }

    *first_filename_idx = optind;

    return kv;
}

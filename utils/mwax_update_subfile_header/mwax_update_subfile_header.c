// 
// Basic utility to replace an existing subfile header with a header based on a supplied text file.
//
// To compile: 
//   gcc -Wall -Ofast mwax_update_subfile_header.c -o mwax_update_subfile_header -lpsrdada
//
// Usage:
//   mwax_update_subfile_header [-h] [-s KEY=VAL [-s ...]] [-d KEY [-d ...]] SUBFILE [SUBFILE ...]
//
//   Run `mwax_update_subfile_header -h` for more detailed usage information.
//
// See the following link for a description of an MWAX Subfile PSRDADA header: https://wiki.mwatelescope.org/display/MP/MWAX+PSRDADA+header
// 
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
#define INCRMT_KV 2 /* Increment value (uses inc_val instead of val) */
struct keyval
{
    char   key[MAX_KV_SIZE];
    char   val[MAX_KV_SIZE];
    int    inc_val;
    int    inc_amount;
    int    instruction; // MODIFY_KV, DELETE_KV, INCRMT_KV
    struct keyval *next;
};

int append_keyval(struct keyval **kv, const char *key, const char *val, int instruction, int inc_amount)
{
    struct keyval *new_kv;

    if (*kv == NULL) // If kv is an empty list...
    {
        *kv = (struct keyval *)malloc(sizeof(struct keyval));
        new_kv = *kv;
    }
    else
    {
        new_kv = *kv;

        // Move to the end of the list
        while (new_kv->next != NULL)
            new_kv = new_kv->next;

        // Create a new node and move to it
        new_kv->next = (struct keyval *)malloc(sizeof(struct keyval));
        new_kv = new_kv->next;
    }

    // Update the values
    if (key != NULL)
        snprintf(new_kv->key, MAX_KV_SIZE, "%s", key);

    if (val != NULL)
        snprintf(new_kv->val, MAX_KV_SIZE, "%s", val);

    new_kv->instruction = instruction;
    new_kv->inc_amount = inc_amount;

    if (instruction == INCRMT_KV)
        new_kv->inc_val = (val == NULL ? 0 : atoi(val));

    new_kv->next = NULL;

    return EXIT_SUCCESS;
}

void usage(FILE *f, char **argv)
{
    fprintf(f, "usage: %s [-h] [-s KEY=VAL [-s ...]] [-d KEY [-d ...]] [-o KEY[=VAL[:INC]]] SUBFILE [SUBFILE ...]\n", argv[0]);
    fprintf(f, "\t-d KEY          Deletes KEY from header.\n");
    fprintf(f, "\t-h              Prints this help and exits.\n");
    fprintf(f, "\t-o KEY[=VAL[:INC]]\n");
    fprintf(f, "\t                Sets KEY=VAL for the first SUBFILE, and increments VAL += INC for\n");
    fprintf(f, "\t                each succesive SUBFILE. If VAL and INC cannot be converted to an int,\n");
    fprintf(f, "\t                or if they are not supplied, they are set to 0 and 1, respectively. INC can be\n");
    fprintf(f, "\t                any integer (including zero and negative). If KEY does not exist, it is created.\n");
    fprintf(f, "\t-s KEY=VAL      Sets the value of VAL to KEY. If KEY does not exist, it is created.\n");
    fprintf(f, "If no -d or -s is given, then the subfile headers are printed out to stdout.\n");
}

struct keyval *parse_cmdline(int argc, char **argv, int *first_filename_idx);

int main(int argc, char **argv)
{
    // Parse the command line
    int first_filename_idx;
    struct keyval *first_kv = parse_cmdline(argc, argv, &first_filename_idx);
    struct keyval *kv;

    // Set up variables for handling I/O
    int file_descr;

    int input_read;
    int input_write;

    char header_buffer[HEADER_LEN];

    // Loop over the remaining arguments, open the files, and view/modify their headers
    int i;
    for (i = first_filename_idx; i < argc; i++)
    {
        printf("Opening %s...\n", argv[i]);

        file_descr = open( argv[i], O_RDWR );
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
        kv = first_kv;
        if (kv == NULL)
        {
            printf("%s:\n%s\n", argv[i], header_buffer);
            close ( file_descr );
            continue;
        }

        // If we got to here, then there are at least some modifications to be made
        while (kv != NULL)
        {
            switch (kv->instruction)
            {
                case MODIFY_KV:
                    if (ascii_header_set(header_buffer, kv->key, "%s", kv->val) == -1)
                    {
                        fprintf(stderr, "WARNING: Could not modify/create %s=%s\n", kv->key, kv->val);
                    }
                    break;
                case DELETE_KV:
                    if (ascii_header_del(header_buffer, kv->key) == -1)
                    {
                        fprintf(stderr, "WARNING: Could not delete keyword %s\n", kv->key);
                    }
                    break;
                case INCRMT_KV:
                    if (ascii_header_set(header_buffer, kv->key, "%d", kv->inc_val) == -1)
                    {
                        fprintf(stderr, "WARNING: Could not increment %s=%d\n", kv->key, kv->inc_val);
                    }
                    kv->inc_val += kv->inc_amount;
                    break;
                default:
                    fprintf(stderr, "WARNING: Unrecognised instruction code (%d)\n", kv->instruction);
            }

            kv = kv->next;
        }

        // Now that the header buffer has been suitably altered, write it back to the beginning of the file
        lseek(file_descr, 0, SEEK_SET);
        input_write = write(file_descr, header_buffer, HEADER_LEN);
        if (input_write != HEADER_LEN)
        {
            printf("FATAL WARNING: Write to %s failed.  Returned a value of %d. "
                    "The header may now be corrupted.\n", argv[i], input_write);
            exit(EXIT_FAILURE);
        }

        close(file_descr);
    }

    exit(EXIT_SUCCESS);
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
    int inc_amount;
    int n; // Used for counting scanf'd items

    int opt;
    while ((opt = getopt(argc, argv, "d:ho:s:")) != -1)
    {
        switch (opt)
        {
            case 'd':
                append_keyval(&kv, optarg, NULL, DELETE_KV, 0);
                break;
            case 'h':
                printf("Modifies subfile headers in-place\n");
                usage(stdout, argv);
                exit(EXIT_SUCCESS);
                break;
            case 'o':
                n = sscanf(optarg, "%[^=]=%[^:]:%d", key, val, &inc_amount);
                if (n == 0)
                {
                    fprintf(stderr, "ERROR: cannot parse \"%s\" as KEY[=VAL[:INC]]\n", optarg);
                    usage(stderr, argv);
                    exit(EXIT_FAILURE);
                }

                if (n == 1)
                    append_keyval(&kv, key, "0", INCRMT_KV, 1);
                else if (n == 2)
                    append_keyval(&kv, key, val, INCRMT_KV, 1);
                else if (n == 3)
                    append_keyval(&kv, key, val, INCRMT_KV, inc_amount);
                else
                {
                    fprintf(stderr, "ERROR: unexpected number (%d) of parsed elements in \"%s\"\n", n, optarg);
                    usage(stderr, argv);
                    exit(EXIT_FAILURE);
                }

                break;
            case 's':
                n = sscanf(optarg, "%[^=]=%s", key, val);
                if (n != 2)
                {
                    fprintf(stderr, "ERROR: cannot parse \"%s\" as KEY=VAL\n", optarg);
                    usage(stderr, argv);
                    exit(EXIT_FAILURE);
                }
                append_keyval(&kv, key, val, MODIFY_KV, 0);
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

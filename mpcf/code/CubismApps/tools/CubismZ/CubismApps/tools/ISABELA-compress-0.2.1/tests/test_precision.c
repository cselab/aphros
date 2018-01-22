#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <float.h>
#include <stdint.h>
#include <math.h>

double get_isabela_double_max_error (double *original, double *approximated, uint32_t n)
{
    uint32_t i;
    double max_error_rate = -DBL_MAX;
    double current_error_rate = 0;

    for (i = 0; i < n; i ++) {

        if (abs (original [i]) <= DBL_EPSILON) continue;
        current_error_rate = fabs ((original [i] - approximated [i]) / original [i] * 100);

        if (current_error_rate > max_error_rate) {
            max_error_rate = current_error_rate; 
        }
    }

    return max_error_rate;
}

float get_isabela_float_max_error (float *original, float *approximated, uint32_t n)
{
    uint32_t i;
    float max_error_rate = -FLT_MAX;
    float current_error_rate = 0;

    for (i = 0; i < n; i ++) {
        if (abs (original [i]) <= FLT_EPSILON) continue;
        current_error_rate = fabs ((original [i] - approximated [i]) / original [i] * 100);

        if (current_error_rate > max_error_rate) {
            max_error_rate = current_error_rate; 
        }
    }

    return max_error_rate;
}

void print_usage_error_exit (const char *_err_msg)
{
    printf ("ERROR: %s!\n", _err_msg);
    printf ("\n");

    printf ("Program Usage Arguments: ");
    printf ("<original file> <decompressed file> <element size>\n");

    printf ("  <original file>  : original data file\n");
    printf ("  <decompressed file> : approximated data file\n");
    printf ("  <element size>: element segment size in bytes (ie, 8 for double, 4 for float)\n");

    exit (EXIT_FAILURE);
}

void read_file (char *file_name, void *input_buffer, size_t input_buffer_size)
{
    FILE* in_file = fopen (file_name, "rb");

    if (in_file == NULL)
        print_usage_error_exit ("Unable to open input file");

    if (fread (input_buffer, 1, input_buffer_size, in_file) != input_buffer_size)
        print_usage_error_exit ("Unable to read input file completely");

    fclose (in_file);
}

int main (int argc, char** argv)
{
    struct stat stat_info;
    size_t original_file_size;
    size_t decompressed_file_size;

    char *original_buffer;
    char *decompressed_buffer;

    int element_bytes;

    if (argc != 4) {
        print_usage_error_exit ("Missing program arguments");
    }

    if (sscanf (argv[3], "%d", &element_bytes) != 1) {
        print_usage_error_exit ("invalid input for element byte size");
    }

    if (stat (argv[1], &stat_info)) {
        printf ("[%s: %d] Unable to read input file %s", __FUNCTION__, __LINE__, argv [1]);
        print_usage_error_exit ("Unable to read original file");
    }

    original_file_size = (size_t) stat_info.st_size;
    original_buffer = (char*) malloc (original_file_size);

    if (stat (argv[2], &stat_info)) {
        printf ("[%s: %d] Unable to read input file %s", __FUNCTION__, __LINE__, argv [1]);
        print_usage_error_exit ("Unable to read decompressed file");
    }
    decompressed_file_size = (size_t) stat_info.st_size;
    decompressed_buffer = (char*) malloc (decompressed_file_size);

    assert (original_buffer != NULL && decompressed_buffer != NULL);

    // if (original_file_size != decompressed_file_size) {
    //     printf ("[%s: %d] Original and decompressed file sizes are not the same.\n", __FUNCTION__, __LINE__);
    //     exit (EXIT_FAILURE);
    // }

    if (original_file_size < decompressed_file_size) {
        decompressed_file_size = original_file_size;
    } else {
        original_file_size = decompressed_file_size;
    }

    uint32_t num_elements = original_file_size / element_bytes;

    // Read original and approximated data
    read_file (argv [1], original_buffer, original_file_size);
    read_file (argv [2], decompressed_buffer, decompressed_file_size);

    // Report the maximum error between the original and the approximated data
    if (element_bytes == sizeof (double)) {
        double max_error = get_isabela_double_max_error ((double *) original_buffer, (double *) decompressed_buffer, num_elements);
        printf ("[%s: %d] Max error between original and approximated data is %lf \%\n", __FUNCTION__, __LINE__, max_error);
    } else if (element_bytes == sizeof (float)) {
        float max_error = get_isabela_float_max_error ((float *) original_buffer, (float *) decompressed_buffer, num_elements);
        printf ("[%s: %d] Max error between original and approximated data is %f \%\n", __FUNCTION__, __LINE__, max_error);
    } else {
        printf ("[%s: %d] Unsupported element type.\n", __FUNCTION__, __LINE__);
        exit (EXIT_FAILURE);
    }

    free (original_buffer);
    free (decompressed_buffer);

    return EXIT_SUCCESS;
}

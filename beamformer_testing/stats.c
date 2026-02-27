#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_ITER 10    // max clipping iterations
#define CLIP_SIGMA 3.0 // clip anything > 3 sigma from mean

// Calculate mean and std of bins marked as off-pulse (mask=1)
static void calc_stats(const double *profile, const int *mask, int nbins,
                       double *mean, double *std)
{
    int n = 0;
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < nbins; i++)
    {
        if (mask[i])
        {
            sum += profile[i];
            sum2 += profile[i] * profile[i];
            n++;
        }
    }
    *mean = sum / n;
    *std = sqrt(sum2 / n - (*mean) * (*mean));
}

double calc_pulsar_snr(const double *profile, int nbins)
{
    printf("DEBUG: nbins = %d\n", nbins);
    printf("DEBUG: first few profile values: %f %f %f %f\n", profile[0], profile[1], profile[2], profile[3]);
    // Start by assuming all bins are off-pulse
    int *mask = malloc(nbins * sizeof(int));
    if (!mask)
        return -1.0;
    for (int i = 0; i < nbins; i++)
        mask[i] = 1;

    // debug - verify mask
    int check = 0;
    for (int i = 0; i < nbins; i++)
        check += mask[i];
    printf("DEBUG: initial mask sum = %d (should be %d)\n", check, nbins);

    double mean, std;

    // Iteratively clip bins that are > CLIP_SIGMA above the mean
    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        calc_stats(profile, mask, nbins, &mean, &std);

        int nclipped = 0;
        for (int i = 0; i < nbins; i++)
        {
            if (mask[i] && (profile[i] - mean) > CLIP_SIGMA * std)
            {
                mask[i] = 0; // mask out this bin (it's on-pulse)
                nclipped++;
            }
        }

        // Converged when no new bins are clipped
        if (nclipped == 0)
            break;
    }

    // Count surviving off-pulse bins
    int noff = 0;
    for (int i = 0; i < nbins; i++)
        noff += mask[i];

    // Subtract baseline from profile and find peak
    double peak = -1e99;
    for (int i = 0; i < nbins; i++)
    {
        double val = profile[i] - mean;
        if (val > peak)
            peak = val;
    }

    double snr = peak / std;

    printf("INFO: off-pulse bins = %d / %d (%.1f%%)\n",
           noff, nbins, 100.0 * noff / nbins);
    printf("INFO: off-pulse mean = %.4f\n", mean);
    printf("INFO: off-pulse std  = %.4f\n", std);
    printf("INFO: peak           = %.4f\n", peak + mean);
    printf("INFO: SNR            = %.2f\n", snr);

    free(mask);

    return snr;
}
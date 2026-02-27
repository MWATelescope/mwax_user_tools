#include <stdio.h>
#include <stdlib.h>

int line_chart_png(const char *filename,
                   const char *title,
                   const char *x_label,
                   const char *y_label,
                   const double *x,
                   const double *y,
                   int n)
{
    FILE *gp = popen("gnuplot", "w");
    if (!gp)
    {
        fprintf(stderr, "ERROR: could not open gnuplot via popen\n");
        return -1;
    }

    // --- Output settings ---
    fprintf(gp, "set terminal pngcairo enhanced font 'Sans,12' size 800,600\n");
    fprintf(gp, "set output '%s'\n", filename);

    // --- Labels and title ---
    fprintf(gp, "set title  '%s'\n", title);
    fprintf(gp, "set xlabel '%s'\n", x_label);
    fprintf(gp, "set ylabel '%s'\n", y_label);

    // --- Style ---
    fprintf(gp, "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 0.5\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key off\n");

    // --- Send data inline via gnuplot's special '-' filename ---
    fprintf(gp, "plot '-' with linespoints ls 1\n");
    for (int i = 0; i < n; i++)
    {
        fprintf(gp, "%f %f\n", x[i], y[i]);
    }
    fprintf(gp, "e\n"); // end of inline data

    pclose(gp);
    return 0;
}

// --- Example usage ---
int example(void)
{
    double x[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    double y[] = {0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100};
    int n = sizeof(x) / sizeof(x[0]);

    if (line_chart_png("chart.png", "y = x^2", "X Values", "Y Values", x, y, n) == 0)
        printf("INFO: chart written to chart.png\n");

    return 0;
}

#include <stdio.h>
#include <omp.h>

const double minReal = -2.0;
const double maxReal = 1.0;
const double minImaginary = -1.0;
const double maxImaginary = 1.0;
const int maxIteration = 1000;
const int m_widht = 550, m_height = 550;

double ConvertReal(int x , int width)
{
    double Real = minReal + (x * (maxReal - minReal) / width);
    return Real;
}

double ConvertImaginary(int y, int height)
{
    double imaginary = minImaginary + (y * (maxImaginary - minImaginary) / height);
    return imaginary;
}

int Mandelbrot(double real, double imag, int Iteration_max)
{
    double zr = 0.0f, zi = 0.0f;
    int count = 0;
    for (count = 0; count < Iteration_max && zr * zr + zi * zi < 4.0; count++)
    {
        double temp = zr * zr - zi * zi + real;
        zi = 2.0 * zr * zi + imag;
        zr = temp;
    }

    return count;
}

int main()
{
    FILE *mandelbrot_img = fopen("MandelbrotOpenMP.ppm", "w");
    fprintf(mandelbrot_img, "P3\n");
    fprintf(mandelbrot_img, "%d %d\n", m_widht, m_height);
    fprintf(mandelbrot_img, "255\n");
    int color[m_height][m_widht];
    double start, end;

    start = omp_get_wtime();

    #pragma omp parallel for
    for (int y = 0; y < m_height; y++)
    {
        #pragma omp parallel for
        for (int x = 0; x < m_widht; x++)
        {
            double real = ConvertReal(x, m_widht);
            double imag = ConvertImaginary(y, m_height);

            color[y][x] = Mandelbrot(real, imag, maxIteration);
        }
    }

    for (int y = 0; y <  m_height; y++)
    {
        for (int x = 0; x < m_widht; x++)
        {
            int red = (color[y][x] % 155);
            int green = (color[y][x] % 128);
            int blue = (color[y][x] % 255);

            fprintf(mandelbrot_img, "%d %d %d\n", red, green, blue);
        }
    }
    end = omp_get_wtime();
    printf("Mandelbrot OMP generado en: %f segundos\n", (end-start));

    fflush(mandelbrot_img);
    fclose(mandelbrot_img);
    return 0;
}

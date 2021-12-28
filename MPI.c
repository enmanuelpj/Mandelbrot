
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include "mpi.h"

const double minReal = -2.0;
const double maxReal = 1.0;
const double minImaginary = -1.0;
const double maxImaginary = 1.0;
const int maxIteration = 1000;
const int m_widht = 550, m_height = 550;

double ConvertReal(int x, int width)
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

int main(int argc, char** argv)
{
    FILE* mandelbrot_img = fopen("MandelbrotMPI.ppm", "w");
    fprintf(mandelbrot_img, "P3\n");
    fprintf(mandelbrot_img, "%d %d\n", m_widht, m_height);
    fprintf(mandelbrot_img, "255\n");
    int color[m_widht];
    int myrank, numProcess, row;
    int Count = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcess);
    MPI_Status status, status2, status3, status4;
    double start, end;
    double syncStart, syncEnd;
    int delta = (m_height / (numProcess - 1)) + 1;

    if (myrank == 0)
    {
        start = MPI_Wtime();

        syncStart = MPI_Wtime();
        for (int i = 0, row = 0; i < numProcess - 1; i++, row = row + delta)
        {
            MPI_Send(&row, 1, MPI_INT, i + 1, i + 1, MPI_COMM_WORLD);
            MPI_Send(&Count, 1, MPI_INT, i + 1, i + 1, MPI_COMM_WORLD);
            MPI_Send(&delta, 1, MPI_INT, i + 1, i + 1, MPI_COMM_WORLD);
            Count = Count + delta;
        }
        syncEnd = MPI_Wtime();
        for (int i = 0; i < m_height; i++)
        {
            MPI_Recv(color, m_widht, MPI_INT, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &status2);
            for (int x = 0; x < m_widht; x++)
            {
                int red = (color[x] % 128);
                int green = (color[x] % 205);
                int blue = (color[x] % 250);
                fprintf(mandelbrot_img, "%d %d %d\n", red, green, blue);
            }
        }

        end = MPI_Wtime();
        printf("Mandelbrot MPI generado en: %f segundos\n", (end - start));
    }
    else
    {
        MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&Count, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status3);
        MPI_Recv(&delta, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status4);
        for (int y = Count; y < row + delta; y++)
        {
            for (int x = 0; x < m_widht; x++)
            {
                double real = ConvertReal(x, m_widht);
                double imag = ConvertImaginary(y, m_height);
                color[x] = Mandelbrot(real, imag, maxIteration);
            }
            MPI_Send(color, m_widht, MPI_INT, 0, y, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
    fflush(mandelbrot_img);
    fclose(mandelbrot_img);
    return 0;
}



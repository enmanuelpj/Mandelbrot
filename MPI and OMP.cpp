#define _CRT_SECURE_NO_DEPRECATE
#include "mpi.h" 
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

#define TAG_ENVIO 3

int ANCHO = 400;
int ALTURA = 400;
int ITERACION = 500;
int PROCESS = 1;
int SIZE;

// Funcion que calcula el valor de la iteracion
void CreacionImagen(char* pixels) {
    int i;
    FILE* fp;
    fp = fopen("Image.ppm", "w"); // Abrimos el archivo

    if (fp == NULL) {
        perror("Error en imagen"); // Error en la imagen
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "P6\n# CREATOR: mandel program\n"); // Escribimos en el archivo
    fprintf(fp, "%d %d\n255\n", ANCHO, ALTURA); // Escribimos en el archivo

    for (i = 0; i < (ANCHO * 3 * ALTURA); i++) {
        fputc((char)pixels[i], fp); // Escribimos en el archivo
    }
    fclose(fp);
}

// Funcion que calcula el valor de la iteracion
void CalcularMandelbrot(char* destino, int iproc, int altura)
{
    // Variables
    int posInici = iproc * altura;
    int xactual, yactual, posLocal = 0;
    int ProcessSize = altura * ANCHO * 3;

    // Calculamos el valor de la iteracion
    omp_set_dynamic(1);
    double pr, pi, nuevoD, nuevoI, viejoD, viejoI;
    double zoom = 1, moveX = -0.5, moveY = 0;
    int numcpu = omp_get_num_procs();
    omp_set_num_threads(numcpu);

    if (iproc != 0) { // Si no es el proceso 0
        destino = (char*)malloc(sizeof(char) * ProcessSize);
    }
    // Calculamos el valor de la iteracion
    #pragma omp parallel 
    {
        #pragma omp for schedule(dynamic) // Para cada iteracion 
        for (yactual = posInici; yactual < posInici + altura; yactual++)
        {
            for (xactual = 0; xactual < ANCHO; xactual++)
            {
                pr = 1.5 * (xactual - ANCHO / 2) / (0.5 * zoom * ANCHO) + moveX;
                pi = (yactual - ALTURA / 2) / (0.5 * zoom * ALTURA) + moveY;
                nuevoD = nuevoI = viejoD = viejoI = 0;
                int i;
                for (i = 0; i < ITERACION; i++)
                {
                    viejoD = nuevoD;
                    viejoI = nuevoI;
                    nuevoD = viejoD * viejoD - viejoI * viejoI + pr;
                    nuevoI = 2 * viejoD * viejoI + pi;
                    if ((nuevoD * nuevoD + nuevoI * nuevoI) > 4) break;
                }

                if (i == ITERACION)
                {
                    destino[posLocal] = 0;
                    destino[++posLocal] = 0;
                    destino[++posLocal] = 0;
                    ++posLocal;
                }
                else
                {
                    double z = sqrt(nuevoD * nuevoD + nuevoI * nuevoI);
                    int brightness = 256 * log2(1.75 + i - log2(log2(z))) / log2((double)ITERACION);
                    destino[posLocal] = brightness;
                    destino[++posLocal] = brightness;
                    destino[++posLocal] = 255;
                    ++posLocal;
                }
            }
        }
    }

    #pragma omp barrier
    if (iproc != 0)
    {
        MPI_Send(destino, ProcessSize, MPI_CHAR, 0, TAG_ENVIO, MPI_COMM_WORLD);
        free(destino);
    }
}

// Funcion que calcula el valor de la iteracion
void stringCopy(char* pixels, char* reciboArr, int sender, int altura)
{
    int posInici = sender * altura * ANCHO * 3;
    int pos;

    for (pos = 0; pos < altura * ANCHO * 3; pos++, posInici++) {
        pixels[posInici] = reciboArr[pos];
    }
}


int main(int argc, char** argv) {
    // Variables
    char* pixels;
    int nproc, iproc;
    double begin, end_total;

    if (argc >= 3) {
        ANCHO = atoi(argv[1]);
        ALTURA = atoi(argv[2]);
    }

    if (argc == 4) {
        ITERACION = atoi(argv[3]);
    }

    if (argc == 5) {
        PROCESS = atoi(argv[4]);
    }

    SIZE = ALTURA * ANCHO * 3; // Tamaño de la imagen

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "Error al inicializar MPI.\n"); // Error al inicializar MPI
        return 100;
    }

    begin = MPI_Wtime();

    if (MPI_Comm_size(MPI_COMM_WORLD, &nproc) != MPI_SUCCESS) {
        fprintf(stderr, "No se puede obtener el contador de procesos.\n");
        MPI_Finalize();
        return 101;
    }
    else if (nproc < 2) {
        fprintf(stderr, "Se necesitan almenos 2 procesos.", nproc);
        MPI_Finalize();
        return 102;
    }

    if (MPI_Comm_rank(MPI_COMM_WORLD, &iproc) != MPI_SUCCESS) {
        fprintf(stderr, "No se puede obtener el rango para el proceso.\n");
        MPI_Finalize();
        return 103;
    }

    if ((ALTURA % nproc) != 0)
    {
        printf("Incompatable numero de procesos.");
        exit(EXIT_FAILURE);
    }

    int altura = ALTURA / nproc;
    int ProcessSize = altura * ANCHO * 3;

    // Inicializamos la imagen
    if (iproc != 0) {
        char* envioArr = (char*)malloc(sizeof(char) * ProcessSize);
        CalcularMandelbrot(envioArr, iproc, altura);
    }
    else if (iproc == 0) {
        pixels = (char*)malloc(SIZE * sizeof(char));
        CalcularMandelbrot(pixels, iproc, altura);
        char* reciboArr = (char*)malloc(sizeof(char) * ProcessSize);
        int i;
        MPI_Status status;

        for (i = 1; i < nproc; i++)
        {
            MPI_Recv(reciboArr, altura * ANCHO * 3, MPI_CHAR, MPI_ANY_SOURCE, TAG_ENVIO, MPI_COMM_WORLD, &status);
            stringCopy(pixels, reciboArr, status.MPI_SOURCE, altura);
        }

        // Guardamos la imagen
        free(reciboArr);
        CreacionImagen(pixels);
        end_total = MPI_Wtime() - begin;
        printf("Tiempo en total: %.10lf segundos \n", end_total); // Tiempo total
        free(pixels);
    }

    MPI_Finalize(); // Finalizamos el MPI
    return 0;
}
using System;
using System.Collections.Generic;
using System.Text;

namespace Algorytmy4
{
    public class GaussP
    {
        private static void StworzMacierz(double[][] A, double[] b)
        {
            for (int k = 0; k < A.Length; k++)
            {
                int numerRzedu = k;
                int numerRzeduMax = ZnajdzRzad(A, k);

                if (numerRzeduMax != numerRzedu)
                {
                    ZamienRzedy(numerRzedu, numerRzeduMax, A, b);
                }
                StworzRzad(A, b, k);
            }
        }

        private static double[] LiczX(double[][] A, double[] b)
        {
            double[] xv = new double[A.Length];
            double pom = 0.0;
            for (int i = b.Length - 1; i >= 0; i--)
            {
                pom = 0.0;
                for (int x = i + 1; x < b.Length; x++)
                {
                    pom = pom + A[i][x] * xv[x];
                }
                xv[i] = (b[i] - pom) / A[i][i];
            }
            return xv;
        }

        private static int ZnajdzRzad(double[][] A, int kolumna)
        {
            int rzad = kolumna;
            int pierwszyrzad = kolumna + 1;
            for (int i = pierwszyrzad; i < A.Length; i++)
            {
                if (A[rzad][kolumna] < A[i][kolumna])
                {
                    rzad = i;
                }
            }
            return rzad;
        }

        private static void ZamienRzedy(int rzad, int rzadmax, double[][] A, double[] b)
        {
            double temp;
            for (int i = 0; i < A.Length; i++)
            {
                temp = A[rzad][i];
                A[rzad][i] = A[rzadmax][i];
                A[rzadmax][i] = temp;
            }

            temp = b[rzad];
            b[rzad] = b[rzadmax];
            b[rzadmax] = temp;
        }

        private static void StworzRzad(double[][] A, double[] b, int k)
        {
            double pom;
            for (int i = k + 1; i < A.Length; i++)
            {
                pom = A[i][k] / A[k][k];
                for (int j = 0; j < A.Length; j++)
                {
                    A[i][j] = A[i][j] - (A[k][j] * pom);


                }
                b[i] = b[i] - (b[k] * pom);
            }
        }
        public static double[] Gauss(double[][] A, double[] bVector)
        {
            StworzMacierz(A, bVector);
            return LiczX(A, bVector);
        }

    }
}

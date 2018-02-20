using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Algorytmy4
{
    public class Aproksymacja
    {
        private int[] rozmiar;
        private double[] wyniki;
        private int koncowy_element;
        private int suma_szeregu;
        private int szereg_sum;

        public Aproksymacja(int[] rozmiar, double[] wyniki, int stopienWielomianu)
        {
            this.rozmiar = rozmiar;
            this.wyniki = wyniki;
            this.suma_szeregu = 2 * stopienWielomianu + 1;
            this.szereg_sum = stopienWielomianu + 1;
            this.koncowy_element = wyniki.Length;
        }

        public double[][] Macierz()
        {
            double[][] matrix = new double[szereg_sum][];

            for (int i = 0; i < matrix.Length; i++)
                matrix[i] = new double[szereg_sum];

            double[] sVector = new double[suma_szeregu];

            for (int x = 0; x < koncowy_element; x++)
                for (int i = 0; i < suma_szeregu; i++)
                    sVector[i] += Math.Pow(rozmiar[x], i);


            for (int x = 0; x < szereg_sum; x++)
                for (int i = 0; i < szereg_sum; i++)
                    matrix[x][i] = sVector[x + i];
            return matrix;
        }

        public double[] Vektor()
        {
            double[] vector = new double[szereg_sum];

            for (int x = 0; x < koncowy_element; x++)
                for (int i = 0; i < szereg_sum; i++)
                    vector[i] += Math.Pow(rozmiar[x], i) * Convert.ToDouble(wyniki.GetValue(x));

            return vector;
        }

    }
}

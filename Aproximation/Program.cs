using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Algorytmy4
{
    class Program
    {

        public static int StopienWielomianu(String nazwa)
        {
            if (nazwa.Equals("gauss"))
                return 3;
            else if (nazwa.Equals("gauss_opt") || nazwa.Equals("seidel"))
                return 2;
            else if (nazwa.Equals("lu"))
                return 1;
            else
                return 0;
        }


        static void Main(string[] args)
        {
            int iloscRownan = 13;
            int iloscKolumn = 7;
            String[] line = new String[iloscKolumn];
            double[][] wyniki = new double[iloscKolumn][];//0-GaussPartial, 1-GaussPartialO, 2-GaussSiedel, 3-Jacobi, 4-EigenPartial, 5-SparseLU
            for (int i = 0; i < wyniki.Length; i++)
            {
                wyniki[i] = new double[iloscRownan];
            }
            try
            {   
                // Open the text file using a stream reader.
            using (StreamReader sr = new StreamReader("C:\\Users\\Billy\\source\\repos\\Algorytmy4\\Algorytmy4\\bezgrzybow.txt"))
            {
                    // Read the stream to a string, and write the string to the console.
                    line = sr.ReadLine().Split(';');
                    double[][] linie = new double[13][];
                    for (int i = 0; i < iloscRownan; i++)
                    {
                        linie[i] = Array.ConvertAll(sr.ReadLine().Split(';'), Double.Parse);
                    }
                
                for(int i =0; i < iloscKolumn; i++)
                    {
                        for (int j = 0; j < iloscRownan; j++)
                        {
                            wyniki[i][j] = linie[j][i];
                        }
                    }
                }
        }
        catch (Exception e)
        {
            Console.WriteLine("The file could not be read:");
            Console.WriteLine(e.Message);
        }
            for (int i = 0; i < iloscKolumn; i++)
            {
                for (int j = 0; j < iloscRownan; j++)
                {
                    Console.WriteLine(wyniki[i][j]);
                }
                Console.WriteLine("---------");
            }
            for (int i = 3; i < iloscKolumn; i++) // główne działanie
            {
                int stopienWielomianu = StopienWielomianu(line[i]); // 3,4,5
                int[] rozmiar = new int[iloscRownan];
                for(int j = 0; j < iloscRownan; j++)
                {
                    rozmiar[j] = Convert.ToInt32(wyniki[1][j]);
                }
                Aproksymacja aproksymacja = new Aproksymacja(rozmiar, wyniki[i], stopienWielomianu);
                double[][] macierz = aproksymacja.Macierz();
                double[] vektor = aproksymacja.Vektor();
                double[] wynik = GaussP.Gauss(macierz, vektor);
                Console.WriteLine("wyniki");
                for (int j = 0; j < wynik.Length; j++)
                {
                    Console.WriteLine(wynik[j]);
                }
            }


            if (System.Diagnostics.Debugger.IsAttached) Console.ReadLine();
        }
    }
}

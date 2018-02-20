#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
//silnia

 int factorial(int n)
 {
    if(n==0) return 1;
    else return n*factorial(n-1);
 }

//wlasna f na potege

double myPow(double x, int n)
 {
  double result = 1;
  int i;
    for (i=0;i<n;i++)
    {
     result = (double)result * x;
    }
    return result;
 }
 
//SZEREG TAYLORA PRZOD

 double sinus(double rad)
 {

   double result = 0;
   int i=0;
    for(i;i<=15;i++)
     {
      result = result + (myPow(-1,i)/factorial(2*i+1))*myPow(rad,(2*i)+1);
     }
    return result;
 }

//SZEREG TAYLORA REVERSE

  double sinusr(double rad)
 {

   double result = 0;
   int i=15;
    for(i;i>=0;i--)
     {
      result = result + (myPow(-1,i)/factorial(2*i+1))*myPow(rad,(2*i)+1);
     }
    return result;
 }

//na podstawie poprzedniego wyrazu(wzor i podstawiam n+1/n - wychodzi mnoznik)

  double sinT(double rad)
  {
    double term = rad;
    double result = rad;
    int i=0;
    for(i;i<10;i++)
    {
     term = term * (-(rad*rad)/((2*i+2)*(2*i+3)));
     result += term;
    }
    return result;
  }


 int main()
 {

    int sum=0, sum2=0, sum3=0;
    double result,i=0;
    double radian;
    
printf("Porownanie dla sin(lib) i liczenie szeregu taylora dla funkcji sinus od przodu: \n");

        for(i;i<90;i++)
        {
            radian = i * M_PI / 180;
            if(sin(radian)==sinus(radian))
            {
                sum++;
            }
            printf("Dla rad: %.20lf \n", radian);
            printf("Moj sin: %.20lf \nStandard sin %.20lf \n\n", sinus(radian),sin(radian));
            printf("Blad bezwzgledny: %.20lf \n", fabs((sinus(radian) - sin(radian))/sin(radian))); 
        }
      
printf("Porownanie dla sin(lib) i liczenie szeregu taylora dla funkcji sinus od tylu: \n");

        for(i=0;i<90;i++)
        {
            radian = i * M_PI / 180;
            if(sin(radian)==sinusr(radian))
            {
                sum2++;
            }
            printf("Dla rad: %.20lf \n", radian);
            printf("Moj sin %.20lf \nStandard sin %.20lf \n\n", sinusr(radian),sin(radian));
            printf("Blad bezwzgledny: %.20lf \n", fabs((sinusr(radian) - sin(radian))/sin(radian))); 
        }
       
printf("Porownanie dla sin(lib) i liczenie wyrazu na podstawie poprzedniego: \n");

        for(i=0;i<90;i++)
        {
            radian = i * M_PI / 180;
            if(sin(radian)==sinT(radian))
            {
                sum3++;
            }
            printf("Dla rad: %.20lf \n", radian);
            printf("Moj sin %.20lf \nStandard sin %.20lf \n\n", sinT(radian),sin(radian));
            printf("Blad bezwzgledny: %.20lf \n", fabs((sinT(radian) - sin(radian))/sin(radian))); 
        }

    printf("Laczna liczba powtorzen A1: %d\n", sum);
    printf("Laczna liczba powtorzen A2: %d\n", sum2);
    printf("Laczna liczba powtorzen A3: %d\n", sum3);
    
   
        
  return 0;
 }

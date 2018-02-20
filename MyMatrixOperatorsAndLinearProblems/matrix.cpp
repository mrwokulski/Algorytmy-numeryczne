//Wiktor Przyby³owski	238230 gr1 inf3

#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <Eigen/Dense>
#include <iomanip>
#include <windows.h>
#define MIN 0.00000000001
#define MAX 99.99999999999
#define MINI 1
#define MAXI 9
using namespace Eigen;
using namespace std;

double PCFreq = 0.0;
__int64 CounterStart = 0;

class Fraction {
    private:
        // Calculates the greates common divisor with
        // Euclid's algorithm
        // both arguments have to be positive
        long long gcd(long long a, long long b) {
            while (a != b) {
                if (a > b) {
                    a -= b;
                } else {
                    b -= a;
                }
            }
            return a;
        }
 
    public:
        long long numerator, denominator;
 
        Fraction() {
            numerator = 0;
            denominator = 1;
        }
 
        Fraction(long long n, long long d) {
            if (d==0) {
                cerr << "Denominator may not be 0." << endl;
                exit(0);
            } else if (n == 0) {
                numerator = 0;
                denominator = 1;
            } else {
                int sign = 1;
                if (n < 0) {
                    sign *= -1;
                    n *= -1;
                }
                if (d < 0) {
                    sign *= -1;
                    d *= -1;
                }
 
                long long tmp = gcd(n, d);
                numerator = n/tmp*sign;
                denominator = d/tmp;
            }
        }
 
        operator int() {return (numerator)/denominator;}
        operator float() {return ((float)numerator)/denominator;}
        operator double() {return ((double)numerator)/denominator;}
};
 
Fraction operator+(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                +rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    return tmp;
}
 
Fraction operator+=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                +rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    lhs = tmp;
    return lhs;
}
 
Fraction operator-(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                -rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    return tmp;
}
 
Fraction operator-=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator
                -rhs.numerator*lhs.denominator,
                lhs.denominator*rhs.denominator);
    lhs = tmp;
    return lhs;
}
 
Fraction operator*(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.numerator,
               lhs.denominator*rhs.denominator);
    return tmp;
}
 
Fraction operator*=(Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.numerator,
               lhs.denominator*rhs.denominator);
    lhs = tmp;
    return lhs;
}
 
Fraction operator*(int lhs, const Fraction& rhs) {
    Fraction tmp(lhs*rhs.numerator,rhs.denominator);
    return tmp;
}
 
Fraction operator*(const Fraction& rhs, int lhs) {
    Fraction tmp(lhs*rhs.numerator,rhs.denominator);
    return tmp;
}
 
Fraction operator/(const Fraction& lhs, const Fraction& rhs) {
    Fraction tmp(lhs.numerator*rhs.denominator,
                 lhs.denominator*rhs.numerator);
    return tmp;
}
 
std::ostream& operator<<(std::ostream &strm, const Fraction &a) {
 
    if (a.denominator == 1) {
        strm << a.numerator;
    } else {
        strm << a.numerator << "/" << a.denominator;
    }
    return strm;
}

Fraction s(0,1);	

template <typename T> class MyMatrix {
private:
    vector<vector<T> > matrix;
    unsigned rows;
    unsigned cols;

    T abs(T x) {
        return (x >= 0) ? x : -x;
    }

    public:
        // konstruktory
        MyMatrix(unsigned wiersz, unsigned kolumna, const T& x) {
            matrix.resize(wiersz);
            for (unsigned i=0; i<matrix.size(); i++) {
                matrix[i].resize(kolumna, x);
            }
            rows = wiersz;
            cols = kolumna;
        }
        MyMatrix(const MyMatrix<T>& pom) {
            matrix = pom.matrix;
            rows = pom.getRowCount();
            cols = pom.getColCount();
        }

        // destruktor
        virtual ~MyMatrix() {}

        // przypisanie
        MyMatrix<T>& operator=(const MyMatrix<T>& pom) {
            if (&pom == this) return *this;
            unsigned new_rows = pom.getRowCount();
            unsigned new_cols = pom.getColCount();

            matrix.resize(new_rows);
            for (unsigned i=0; i<matrix.size(); i++) {
                matrix[i].resize(new_cols);
            }

            for (unsigned i=0; i<new_rows; i++) {
                for (unsigned j=0; j<new_cols; j++) {
                    matrix[i][j] = pom(i, j);
                }
            }
            rows = new_rows;
            cols = new_cols;
            return *this;
        }

        // dodawanie macierzy
        MyMatrix<T> operator+(const MyMatrix<T>& pom) {
            MyMatrix result(rows, cols, s);
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    result(i,j) = this->matrix[i][j] + pom(i,j);
                }
            }
            return result;
        }
        MyMatrix<T>& operator+=(const MyMatrix<T>& pom) {
            unsigned rows = pom.getRowCount();
            unsigned cols = pom.getColCount();
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    this->matrix[i][j] += pom(i,j);
                }
            }
            return *this;
        }
        MyMatrix<T> operator-(const MyMatrix<T>& pom) {
            MyMatrix result(rows, cols, s);
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    result(i,j) = this->matrix[i][j] - pom(i,j);
                }
            }
            return result;
        }
        MyMatrix<T>& operator-=(const MyMatrix<T>& pom) {
            unsigned rows = pom.getRowCount();
            unsigned cols = pom.getColCount();
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    this->matrix[i][j] -= pom(i,j);
                }
            }
            return *this;
        }
        MyMatrix<T> operator*(const MyMatrix<T>& pom) {
            unsigned rows = pom.getRowCount();
            unsigned cols = pom.getColCount();
            MyMatrix result(rows, cols, s);
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    for (unsigned k=0; k<rows; k++) {
                        result(i,j) += this->matrix[i][k] * pom(k,j);
                    }
                }
            }
            return result;
        }
        MyMatrix<T>& operator*=(const MyMatrix<T>& pom) {
            MyMatrix result = (*this) * pom;
            (*this) = result;
            return *this;
        }
        MyMatrix<T> transpose() {
            MyMatrix result(rows, cols, s);
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    result(i,j) = this->matrix[j][i];
                }
            }
            return result;
        }

        // skalary
        MyMatrix<T> operator+(const T& pom) {
            MyMatrix result(rows, cols, s);
            for (unsigned i = 0; i < rows; i++) {
                for (unsigned j = 0; j < cols; j++) {
                    result(i,j) = this->matrix[i][j] + pom;
                }
            }
            return result;
        }
        MyMatrix<T> operator-(const T& pom) {
            MyMatrix result(rows, cols, s);
            for (unsigned i = 0; i < rows; i++) {
                for (unsigned j = 0; j < cols; j++) {
                    result(i,j) = this->matrix[i][j] - pom;
                }
            }
            return result;
        }
        MyMatrix<T> operator*(const T& pom) {
            MyMatrix result(rows, cols, s);
            for (unsigned i = 0; i < rows; i++) {
                for (unsigned j = 0; j < cols; j++) {
                    result(i,j) = this->matrix[i][j] * pom;
                }
            }
            return result;
        }
        MyMatrix<T> operator/(const T& pom) {
            MyMatrix result(rows, cols, s);
            for (unsigned i = 0; i < rows; i++) {
                for (unsigned j = 0; j < cols; j++) {
                    result(i,j) = this->matrix[i][j] / pom;
                }
            }
            return result;
        }

        // mnozenie macierzy przez wektor
        vector<T> operator*(const vector<T>& pom) {
            vector<T> result(pom.size(), s);
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    result[i] += this->matrix[i][j] * pom[j];
                }
            }
            return result;
        }
        vector<T> diagonalVector() {
            vector<T> result(rows, s);
            for (unsigned i=0; i<rows; i++) {
                result[i] = this->matrix[i][i];
            }
            return result;
        }

        // dostan sie do poszczególnych elementóœ
        T& operator()(const unsigned& row, const unsigned& col) {
            return this->matrix[row][col];
        }
        const T& operator()(const unsigned& row, const unsigned& col) const {
            return this->matrix[row][col];
        }

        // gettery i settery do kolumn
        void setAt(unsigned row, unsigned col, const T& x) {
            this->matrix[row][col] = x;
        }

        T& getAt(unsigned row, unsigned col) {
            return this->matrix[row][col];
        }

        // zwróc liczbe wierszy/kolumn
        unsigned getRowCount() const {
            return this->rows;
        }
        unsigned getColCount() const {
            return this->cols;
        }

        // wyswietl
        void display() {
            for (int i = 0; i < getRowCount(); i++) {
                cout << "[ ";
                for (int j=0; j < getColCount(); j++) {
                    cout << setw(8) << getAt(i,j);
                }
                cout << " ]" << endl;
            }
        }
        // metoda Gauss-Jordan pelnego wyboru
        vector<T> Gauss3() {
            int n = getRowCount();
            for (int i = 0; i < n; i++) {
                // znajdz wiersz z maksymalnym elementem
                T maxEl = abs(matrix[i][i]);
                int maxRow = i;
                int maxCol = i;
                for (int k = i+1; k < n; k++) {
                    if (abs(matrix[k][i]) > maxEl) {
                        maxEl = abs(matrix[k][i]);
                        maxRow = k;
                        maxCol = i;
                    }
                }
                // zamien maksymalny wiersz z obecnym
                for (int k = i; k < n+1; k++) {
                    T pom = matrix[maxRow][k];
                    matrix[maxRow][k] = matrix[i][k];
                    matrix[i][k] = pom;
                }
                // zamien maksymalna kolumne z obecna
                for (int k = i; k < n; k++) {
                    T pom = matrix[k][maxCol];
                    matrix[k][maxCol] = matrix[k][i];
                    matrix[k][i] = pom;
                }
               
                // wyprowadz zera przed obecnym wierszem
                for (int k = i+1; k < n; k++) {
                    T c = -matrix[k][i] / matrix[i][i];
                    for (int j = i; j < n+1; j++) {
                        if (i == j) {
                            matrix[k][j] = 0;
                        } else {
                            matrix[k][j] += c * matrix[i][j];
                        }
                    }
                }
            }
            // rozwiaz Ax = B za pomoca powstalej macierzy trójkatnej
            vector<T> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            return x;
        }
 
 
 
   // metoda Gauss-Jordan polowicznego wyboru (tylko wiersze)
        vector<T> Gauss2() {
            int n = getRowCount();
            for (int i = 0; i < n; i++) {
                // znajdz wiersz z maksymalnym elementem
                T maxEl = abs(matrix[i][i]);
                int maxRow = i;
                for (int k = i+1; k < n; k++) {
                    if (abs(matrix[k][i]) > maxEl) {
                        maxEl = abs(matrix[k][i]);
                        maxRow = k;
                    }
                }
                // zamien maksymalny wiersz z obecnym
                for (int k = i; k < n+1; k++) {
                    T pom = matrix[maxRow][k];
                    matrix[maxRow][k] = matrix[i][k];
                    matrix[i][k] = pom;
                }
                // wyprowadz zera przed obecnym wierszem
                for (int k = i+1; k < n; k++) {
                    T c = -matrix[k][i] / matrix[i][i];
                    for (int j = i; j < n+1; j++) {
                        if (i == j) {
                            matrix[k][j] = 0;
                        } else {
                            matrix[k][j] += c * matrix[i][j];
                        }
                    }
                }
            }
            // rozwiaz Ax = B za pomoca powstalej macierzy trójkatnej
            vector<T> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            return x;
        }
       
        //bez wyboru
          vector<T> Gauss() {
            int n = getRowCount();
            for (int i = 0; i < n; i++) {
                             
               
                // wyprowadz zera przed obecnym wierszem
                for (int k = i+1; k < n; k++) {
                    T c = -matrix[k][i] / matrix[i][i];
                    for (int j = i; j < n+1; j++) {
                        if (i == j) {
                            matrix[k][j] = 0;
                        } else {
                            matrix[k][j] += c * matrix[i][j];
                        }
                    }
                }
            }
        }
};

double dRand(double min, double max)
{
    double f = (double)rand() / RAND_MAX;
    return min + f * (max - min);
}

int Rand(int max)
{
	int min = 1;
   int output;
   output = min + (rand() % static_cast<int>(max - min + 1));
   return output;
}

double D_(Fraction l){
	double wynik;
	wynik = (double)l.numerator/l.denominator;
	return wynik;	
}

void StartCounter()
{
    LARGE_INTEGER li;
    if(!QueryPerformanceFrequency(&li))
    cout << "QueryPerformanceFrequency failed!\n";

    PCFreq = double(li.QuadPart)/1000000.0;

    QueryPerformanceCounter(&li);
    CounterStart = li.QuadPart;
}
double GetCounter()
{
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    return double(li.QuadPart-CounterStart)/PCFreq;
}

double blad(double f1, double f2){
	double wynik = f1 - f2;
	if(wynik < 0) wynik *= (-1);
	return wynik;
}

int main(int argc, char** argv) {

	srand(time(0));
	
	FILE *plik = fopen("ax.csv", "w");
   fprintf(plik, "eigen;double;float\n");
	
	int y;
	cout << "Podaj wielkosc";
	cin >> y;

	for(int i=2; i<=y; i++){
		cout << "\n";
		
		int ilosc = i*i;
		int ulamek[ilosc][2];
		float liczba1[ilosc];
		double liczba2[ilosc];
		float wektor_f[i];
		double wektor_d[i];
		
		MyMatrix<Fraction> A_fract(i,i,s);
		MyMatrix<float> A_float(i,i+1,1.0);
		MyMatrix<double> A_double(i,i+1,1.0);
		
		MyMatrix<Fraction> B_fract(i,i,s);
		MyMatrix<float> B_float(i,i+1,1.0);
		MyMatrix<double> B_double(i,i+1,1.0);
		
		MyMatrix<Fraction> C_fract(i,i,s);
		MyMatrix<float> C_float(i,i,1.0);
		MyMatrix<double> C_double(i,i,1.0);
		
		MatrixXd Eigen_A = MatrixXd(i,i);
		MatrixXd Eigen_B = MatrixXd(i,i);
		MatrixXd Eigen_C = MatrixXd(i,i);
		
		
		vector<Fraction> X_frac(i);
		vector<float> X_float(i);
		vector<double> X_double(i);
		VectorXd Eigen_X(i);
		
		//tworzê u³amki do macierzy	
		Fraction ulamki[ilosc]; 	
		
		for(int j=0; j<(ilosc); j++){
		//	cout << " i = " << ilosc << " j =" << j << endl;		
			ulamek[j][0] = Rand(10);
			ulamek[j][1] = Rand(10);
			liczba1[j] = (float) ulamek[j][0] / ulamek[j][1];
			liczba2[j] = (double) ulamek[j][0] / ulamek[j][1];
		}
		
			for(int j=0; j<(ilosc); j++){					
			Fraction u(ulamek[j][0],ulamek[j][1]);
			ulamki[j] = u;
		}
		//tworzê wektor
		int wektor_ulamki[ilosc][2];
		Fraction wektor_ulamkow[ilosc];
		
			for(int j=0; j<i; j++){
			wektor_ulamki[j][0] = Rand(10);
			wektor_ulamki[j][1] = Rand(10);
			wektor_f[j] = (float) wektor_ulamki[j][0] / wektor_ulamki[j][1];
			wektor_d[j] = (double) wektor_ulamki[j][0] / wektor_ulamki[j][1];			
		}
		
		for(int j=0; j<i; j++){					
			Fraction u2(wektor_ulamki[j][0],wektor_ulamki[j][1]);
			wektor_ulamkow[j] = u2;
		}
		
		//wrzucam do matrixa		  
   			 int p = 0;
   			  int d = 0;
			  for (int j=0; j<i; j++) {
		      for(int c=0; c<i;c++){
			   //cout << p << " " << c << " j =" << d << endl ;
			  A_fract(p,c) = ulamki[d];
			
			  
			  A_double(p,c) = liczba2[d];
			  A_float(p,c) = liczba1[d];
			  Eigen_A(p,c) = liczba2[d];

			  
			  d++;
			  }
			 p++;			 		     
		    }
		    
		  //wrzucam do wektora  
		  for (int j=0; j<i; j++){		  
 		   X_frac[j] = wektor_ulamkow[j];
 		   X_float[j] = wektor_f[j];
 		   X_double[j] = wektor_d[j];
 		   Eigen_X[j] = wektor_d[j];
 		  
 		  }

		    for(int c=0; c<A_double.getRowCount(); c++){
  			A_double(c,A_double.getColCount()-1)=wektor_d[c];
  			}
  			
  			for(int c=0; c<A_float.getRowCount(); c++){
  			A_float(c,A_float.getColCount()-1)=wektor_f[c];
  			}
  				 	
		 	VectorXd wynik = Eigen_A.partialPivLu().solve(Eigen_X);     
		    vector<double> wynik_duble = A_double.Gauss2() ;		     
		    vector<float> wynik_float = A_float.Gauss2();
		
					for(int j=0; j< i; j++){
						cout << "wynikowe " << j << endl;
					fprintf(plik, "%.15lf;", blad(wynik[j],wynik[j]));
					fprintf(plik, "%.15lf;",  blad(wynik[j],wynik_duble[j]));
					fprintf(plik, "%.15lf\n",  blad(wynik[j],wynik_float[j]));
				}
									
				
	}

fclose(plik);
	
 
    
    return 0;
}
  
  
  
  
/*	
	cout.precision(17);
    int x = 3,m,n,tmp;
	srand(time(0));

//for(int i=2 ; i<x; i++){

    MatrixXd eA = MatrixXd(x,x);
    MatrixXd eAgauss = MatrixXd(x, x);
    MatrixXd eB = MatrixXd(x,x);
    MatrixXd eC = MatrixXd(x,x);
    VectorXd evec = VectorXd(x);
//================================================ DOUBLE

    double MatrixA[x][x];
    double MatrixB[x][x];
    double MatrixC[x][x];
    double vec[x];
    
    MyMatrix<double> A(x, x, 1);
    MyMatrix<double> Agauss(x, x+1, 1);
    MyMatrix<double> B(x, x, 1);   
    MyMatrix<double> C(x, x, 1);
    vector<double> vecd(x);
    vector<double> wynd(x);
    
//wypelnianie tablic rand wartosciami
    
	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
    		tmp = dRand(MIN,MAX);
      		MatrixA[i][j] = tmp;
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
    		tmp = dRand(MIN,MAX);
      		MatrixB[i][j] = tmp; 	
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
    		tmp = dRand(MIN,MAX);
      		MatrixC[i][j] = tmp;
    	}  
  	}
  	for(int i=0; i<x; i++){
  		vec[i]= dRand(MIN,MAX);
	}
    
//wstawianie tablicy do ciala MyMatrix i Eigen
    
	for (int i=0; i<x; i++) {
		tmp = dRand(MIN,MAX);
		Agauss(i,x) = tmp;
    	for (int j=0; j<x; j++) {
      		A(i,j) = MatrixA[i][j];
      		Agauss(i,j) = MatrixA[i][j];
      		eA(i,j) = MatrixA[i][j];
      		eAgauss(i,j) = MatrixA[i][j];
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
      		B(i,j) = MatrixB[i][j];
      		eB(i,j) = MatrixB[i][j];
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
      		C(i,j) = MatrixC[i][j];
      		eC(i,j) = MatrixC[i][j];
    	}  
  	}
  	for(int i=0; i<x; i++){
  		vecd[i]=vec[i];
  		evec[i]=vec[i];
	}

cout << "=====================" << endl;
VectorXd ewyn1 = eA*evec;
cout << ewyn1 << endl;

//1	
	StartCounter();
	vector<double>wynd1 = A*vecd;
	cout << GetCounter()<< " ";
	
	for(int i=0 ; i<x ; i++){
	cout << wynd1[i] << endl;
}
	
	
//2	
	StartCounter();
	MyMatrix<double> h = A+B+C;
	vector<double> wynd2 = h*vecd;
	cout << GetCounter()<< " ";	
//3	
	StartCounter();
	MyMatrix<double> wynd3 = A*(B*C);
	cout << GetCounter()<< " ";
//4	
	StartCounter();
	vector<double> wyndg1 = Agauss.Gauss1();
	cout << GetCounter() << " ";
//5
	StartCounter();
	vector<double> wyndg2 = Agauss.Gauss2();
	cout << GetCounter() << " ";
//6	
	StartCounter();
	vector<double> wyndg3 = Agauss.Gauss3();
	cout << GetCounter() << " ";
	StartCounter();	
	eA.fullPivLu().solve(evec);
	cout << GetCounter() << endl;
//1	eigen
	StartCounter();
	VectorXd ewyn1 = eA*evec;
	cout << GetCounter() << endl;
//2 
	StartCounter();
	VectorXd ewyn2 = (eA+eB+eC)*evec;
	cout << GetCounter() << endl;		
//3	
	StartCounter();
	MatrixXd ewyn3 = eA*(eB*eC);
	cout << GetCounter() << endl;
//4 
	StartCounter();	
	eA.partialPivLu().solve(evec);
	cout << GetCounter() << endl;
//5	
	StartCounter();	
	eA.fullPivLu().solve(evec);
	cout << GetCounter() << endl;


//=================================================== FLOAT
/*	
	float MatrixA[x][x];
    float MatrixB[x][x];
    float MatrixC[x][x];
    float vec[x];
    
    MyMatrix<float> A(x, x, 1);
    MyMatrix<float> Agauss(x, x+1, 1);
    MyMatrix<float> B(x, x, 1);   
    MyMatrix<float> C(x, x, 1);
    vector<float> vecd(x);
    vector<float> wynd(x);

//wypelnianie tablic rand wartosciami
    
	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
    		tmp = dRand(MIN,MAX);
      		MatrixA[i][j] = tmp;
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
    		tmp = dRand(MIN,MAX);
      		MatrixB[i][j] = tmp; 	
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
    		tmp = dRand(MIN,MAX);
      		MatrixC[i][j] = tmp;
    	}  
  	}
  	for(int i=0; i<x; i++){
  		vec[i]= dRand(MIN,MAX);
	}
    
//wstawianie tablicy do ciala MyMatrix i Eigen
    
	for (int i=0; i<x; i++) {
		tmp = dRand(MIN,MAX);
		Agauss(i,x) = tmp;
    	for (int j=0; j<x; j++) {
      		A(i,j) = MatrixA[i][j];
      		Agauss(i,j) = MatrixA[i][j];
      		eA(i,j) = MatrixA[i][j];
      		eAgauss(i,j) = MatrixA[i][j];
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
      		B(i,j) = MatrixB[i][j];
      		eB(i,j) = MatrixB[i][j];
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
      		C(i,j) = MatrixC[i][j];
      		eC(i,j) = MatrixC[i][j];
    	}  
  	}
  	for(int i=0; i<x; i++){
  		vecd[i]=vec[i];
  		evec[i]=vec[i];
	}
	
//1	
	StartCounter();
	vector<float>wynd1 = A*vecd;
	cout << GetCounter()<< " ";
	StartCounter();
	VectorXd ewyn1 = eA*evec;
	cout << GetCounter() << endl;
//2	
	StartCounter();
	MyMatrix<float> h = A+B+C;
	vector<float> wynd2 = h*vecd;
	cout << GetCounter()<< " ";	
	StartCounter();
	VectorXd ewyn2 = (eA+eB+eC)*evec;
	cout << GetCounter() << endl;
//3	
	StartCounter();
	MyMatrix<float> wynd3 = A*(B*C);
	cout << GetCounter()<< " ";
	StartCounter();
	MatrixXd ewyn3 = eA*(eB*eC);
	cout << GetCounter() << endl;
//4	
	StartCounter();
	vector<float> wyndg1 = Agauss.Gauss();
	cout << GetCounter() << endl;
//5
	StartCounter();
	vector<float> wyndg2 = Agauss.Gauss2();
	cout << GetCounter() << " ";
	StartCounter();	
	eA.partialPivLu().solve(evec);
	cout << GetCounter() << endl;
//6	
	StartCounter();
	vector<float> wyndg3 = Agauss.Gauss3();
	cout << GetCounter() << " ";
	StartCounter();	
	eA.fullPivLu().solve(evec);
	cout << GetCounter() << endl;
//1	eigen
	StartCounter();
	VectorXd ewyn1 = eA*evec;
	cout << GetCounter() << endl;
//2 
	StartCounter();
	VectorXd ewyn2 = (eA+eB+eC)*evec;
	cout << GetCounter() << endl;		
//3	
	StartCounter();
	MatrixXd ewyn3 = eA*(eB*eC);
	cout << GetCounter() << endl;
//4 
	StartCounter();	
	eA.partialPivLu().solve(evec);
	cout << GetCounter() << endl;
//5	
	StartCounter();	
	eA.fullPivLu().solve(evec);
	cout << GetCounter() << endl;
*/
//===================================================
/*	
 	MyMatrix<Fraction> Af(x, x, s);
    MyMatrix<Fraction> Afgauss(x, x+1, s);
    MyMatrix<Fraction> Bf(x, x, s);   
    MyMatrix<Fraction> Cf(x, x, s);
	vector<Fraction> vecf(x);
	vector<Fraction> wynf(x);
	
	for (int i=0; i<x; i++) {
		m = Rand(MINI,MAXI);
      	n = Rand(MINI,MAXI);
      	Fraction fill(m,n);
		Afgauss(i,x)=fill;
		evec[i]=(double)fill;
    	for (int j=0; j<x; j++) {
      		m = Rand(MINI,MAXI);
      		n = Rand(MINI,MAXI);
      		Fraction fill(m,n);
     		Af(i,j)=fill;
     		Afgauss(i,j)=fill;
     		eA(i,j) = (double)fill;
      		eAgauss(i,j) = (double)fill;
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
      		m = Rand(MINI,MAXI);
      		n = Rand(MINI,MAXI);
      		Fraction fill(m,n);
     		Bf(i,j)=fill;
     		eB(i,j) = (double)fill;
    	}  
  	}
  	for (int i=0; i<x; i++) {
    	for (int j=0; j<x; j++) {
      		m = Rand(MINI,MAXI);
      		n = Rand(MINI,MAXI);
      		Fraction fill(m,n);
     		Cf(i,j)=fill;
     		eC(i,j) = (double)fill;
    	}  
  	}
  	for(int i=0; i<x; i++){
  		m = Rand(MINI,MAXI);
      	n = Rand(MINI,MAXI);
      	Fraction fill(m,n);
  		vecf[i]=fill;
  		evec[i]=(double)fill;
	}

//1	frac

	StartCounter();
	vector<Fraction> wyn1 = Af*vecf;
	cout << GetCounter() << " ";
	
//2
	StartCounter();
	vector<Fraction> wyn2 = (Af + Bf + Cf)*vecf;
	cout << GetCounter() <<  " ";
	
//3
	StartCounter();
	MyMatrix<Fraction> h = Bf*Cf;
 	MyMatrix<Fraction> wyn3 = Af*h;
 	cout << GetCounter() <<  " ";

//1	eigen
	StartCounter();
	VectorXd ewyn1 = eA*evec;
	cout << GetCounter() << endl;
	
//2 
	
	StartCounter();
	VectorXd ewyn2 = (eA+eB+eC)*evec;
	cout << GetCounter() << endl;
		
//3	
	StartCounter();
	MatrixXd ewyn3 = eA*(eB*eC);
	cout << GetCounter() << endl;
//4 
	StartCounter();	
	eA.partialPivLu().solve(evec);
	cout << GetCounter() << endl;
//5	
	StartCounter();	
	eA.fullPivLu().solve(evec);
	cout << GetCounter() << endl;
*/


//}


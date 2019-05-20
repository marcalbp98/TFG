#include <iostream>
using namespace std;
//Defining the adaptation rate ----- adaptation -> a ,  adaptation rate -> d_a
double d_a(double a,double r, double t){
	d_a= (-a) + r;
	return d_a;
}
//Defining the firing rate ------ firing rate of the population -> d_r , number of cells firing -> r
double d_r(double r, double a, double phi ,double t){
	d_r= -r + a*phi;
	return d_r;
}
int main(){
	double r = 0;
	double a = 0;
	double phi = 0;
	double t = 0;

	cout << "introduce the number of cells firing 'r' " << endl;
	cin << r;
	cout << "introduce the cells that adapt in the population" << endl;
	cin << a;
	cout << "introduce the phi" << endl;
	cin << phi;
	cout << "introduce the time" << endl;
	cin << t;
	d_a(a,r,t);
	d_r(r,a,phi,t);
}

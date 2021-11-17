
using namespace std;
#include <array>
#include <iostream>
#include <string>
#define M_PI 3.14
#define N 64
#define r 1.0
#define ro 0.68

#include <cmath>


const double l = r * 2 * sqrt(N);
const double a = ro * N;
const double b = 2 * ro * l * sqrt(N);
const double c = ro * pow(l, 2) - M_PI * pow(r, 2) * N;
const double D = pow(b, 2) - 4 * a * c;
const double ds = (-b + sqrt(D)) / (2 * a); //расстояние между частицами
const double L = l + sqrt(N) * ds; //размер поля
const double S = pow(L, 2) - N * pow(ds, 2); //площадь, занимаемая частицами


class Particle {
public:
	double x;
	double y;
	//nearest particle
	int index;
	//int index_of_nearest;
	array<int, N - 1> x_nearest_indeces; //массив индексов соседей по возрастанию расстояния
	array<int, N - 1> y_nearest_indeces; //массив индексов соседей по возрастанию расстояния
	bool flag; // particle participated in Event-Chain or not
	Particle() {
		x = 0;
		y = 0;
		index = -1;
		flag = 0;
	}
	Particle(double x_, double y_, int index_) {
		x = x_;
		y = y_;
		index = index_;
		flag = 0;
	}
	/*
	Particle* get_near(Particle A[], Particle a) {
		//Particle B(A.x_nearest, A.y_nearest, index_of_nearest);
		//return &B;
		return &A[a.index_of_nearest];
	}
	void record_nearest(Particle B) {
		x_nearest = B.x;
		y_nearest = B.y;
		index_of_nearest = B.index;
	}
		void get_nearest(double x_n, double y_n, int index_) {
		x_n = x_nearest;
		y_n = y_nearest;
		index_ = index_of_nearest;

	}*/
	void record_new_position(double x_, double y_) {
		x = x_;
		y = y_;
	}

	friend void print(Particle A);
	friend void print_all(Particle A[]);


};

void print(Particle A) {
	cout << "(" << A.x << ", " << A.y << ")";
}

void print_all(Particle A[]) {
	for (int i = 0; i < N; i++) {
		print(A[i]);
		if (i != N - 1) cout << ", ";
		else cout << endl;
	}
	cout << "x = [";
	for (int i = 0; i < N; i++) {
		cout << A[i].x;
		if (i != N - 1) cout << ", ";
		else cout << "]" << endl;
	}
	cout << "y = [";
	for (int i = 0; i < N; i++) {
		cout << A[i].y;
		if (i != N - 1) cout << ", ";
		else cout << "]" << endl;
	}
}
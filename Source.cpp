using namespace std;
#include "Particle.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <vector>
#include <string>

#include <fstream>
#include <math.h>
double get_distance_y(const Particle& A, const Particle& B) {
	double dx = abs(B.x - A.x);

	double dy = B.y - A.y;
	if (dy < 0) dy += L;

	 // check if particle is reachable
	if (dx > 2 * r) return numeric_limits<double>::max();
	return dy - sqrt(4 * r * r - dx * dx);
}
double get_distance_x(const Particle& A, const Particle& B) {
	double dy = abs(B.y - A.y);

	double dx = B.x - A.x;
	if (dx < 0) dx += L;
	
	// check if particle is reachable
	if (dy > 2 * r) return numeric_limits<double>::max();
	
	return dx - sqrt(4 * r * r - dy * dy);
}
void update_nearest(array<Particle, N>& A)
{
	for (int i = 0; i != N; i++)
	{
		std::sort(A[i].y_nearest_indeces.begin(), A[i].y_nearest_indeces.end(),
			[i, A](int lhs, int rhs)
			{
				return get_distance_y(A[i], A[lhs]) < get_distance_y(A[i], A[rhs]);
			}
		);
		std::sort(A[i].x_nearest_indeces.begin(), A[i].x_nearest_indeces.end(),
			[i, A](int lhs, int rhs)
			{
				return get_distance_x(A[i], A[lhs]) < get_distance_x(A[i], A[rhs]);
			}
		);
	}
}
void creating_field(array<Particle, N>& p) {
	int floor = 0, index = 0;
	for (double j = r; j < L; j += (2 * r + ds)) {
		for (double i = r; i < L; i += (2 * r + ds)) {
			p[index].x = i + floor % 2 * r;
			p[index].y = j;
			index++;
		}
		floor++;
	}

	// init nearest
	for (int i = 0; i != N; i++)
	{
		Particle* a = &p[i];
		a->index = i;
		bool magic_flag = false;
		for (int j = 0; j != N; j++)
		{
			if (magic_flag)
			{
				a->x_nearest_indeces[j - 1] = j;
				a->y_nearest_indeces[j - 1] = j;
			}
			else
			{
				if (j == i)
				{
					magic_flag = true;
					continue;
				}
				a->x_nearest_indeces[j] = j;
				a->y_nearest_indeces[j] = j;
			}

		}
	}

	// update nearest
	update_nearest(p);
}
void collision(bool dir, Particle* A, double delta) {
	if (dir)
	{
		A->y += delta;
		if (A->y >= L) A->y -= L;
	}
	else
	{
		A->x += delta;
		if (A->x >= L) A->x -= L;
	}
}
bool end_of_onechain(double distance) {
	if (distance < sqrt(N) * 2 * r)
		return 0;
	else
		return 1;
}
bool is_reachable(const Particle& a, const Particle& b, bool dir)
{
	if (dir) return get_distance_y(a, b) != numeric_limits<double>::max();
	else return get_distance_x(a, b) != numeric_limits<double>::max();
}
Particle* get_right_nearest(array<Particle, N>& A, Particle* a, bool direction) {
	if (direction) return &A[a->y_nearest_indeces[0]];
	else return &A[a->x_nearest_indeces[0]];
}
void change_position(Particle* A, double distance, bool direction) { //перемещение частицы А
	//используется, когда оставшееся расстояние для частицы в цикле в OneChain меньше
	//расстояния до ближайшей в заданном направлении.
	if (direction == 0) {
		A->x += distance;
		if (A->x >= L)
			A->x -= L;
	}
	else {
		A->y += distance;
		if (A->y >= L)
			A->y -= L;
	}
}
ofstream out("particles.txt", ios::out); // открываем файл для записи
void recording(array<Particle, N>& A, ofstream& out) {
	if (out.is_open()) {
		for (int i = 0; i < N; i++) {
			out << A[i].x;
			if (i != N - 1) out << " ";
			else out << endl;
		}
		for (int i = 0; i < N; i++) {
			out << A[i].y;
			if (i != N - 1) out << " ";
			else out << endl;
		}
	}
}
void OneChain(array<Particle, N>& A, int index = 0) { //dir - direction: 0 - right, 1 - up
													   //index - index of start particle, =0
	double distance = 0; //сколько всего пройдено
	bool direction;
	double delta = 0.0;
	direction = rand() % 2;
	Particle* a = &A[index];
	Particle* B;
	while (end_of_onechain(distance) == 0) {
		B = get_right_nearest(A, a, direction);
		if (direction) delta = get_distance_y(*a, *B);
		else delta = get_distance_x(*a, *B);

		if (distance + delta > sqrt(N) * 2 * r) {
			change_position(a, sqrt(N) * 2 * r - distance, direction);
			break;
		}
		if (delta < 0.01) break;
		collision(direction, a, delta);
		distance += delta;
		update_nearest(A);
		a = B;
		recording(A, out);
	}

}
vector<double> rad_func(const array<Particle, N>& p)
{
	//histo
	vector<int> hist;
	double dr = 0.01*r;
	for (double d = 0; d <= 0.1; d += dr/r)
	{
		hist.push_back(0);
		for (int i = 0; i < N - 1; i++)
		{
			for (int j = i + 1; j < N; j++)
			{
				/* double dx = p[j].x - p[i].x;
				double dy = p[j].y - p[i].y; */
				if (get_distance_y(p[i], p[j]) * get_distance_y(p[i], p[j]) + get_distance_x(p[i], p[j]) * get_distance_x(p[i], p[j]) <= pow(d * r + dr/ 2, 2))
				{
					if (get_distance_y(p[i], p[j]) * get_distance_y(p[i], p[j]) + get_distance_x(p[i], p[j]) * get_distance_x(p[i], p[j]) >= pow(d * r - dr / 2, 2))
						hist[d * 100]++;
				}
			}
		}
	}

	
	//g
	vector<double> g(11);
	for (int j = 0; j < 11; j++)
		if (hist[j] != 0)
			g[j] = L * L * hist[j] / (N * 2 * M_PI * (2 * r + j * 0.01) * 0.01);

	return g;
}
double avarage(const vector<double>& g)
{
	if (g.empty()) return 0.0;
	double result = 0;
	size_t t = 1;
	for (auto g_i : g)
		result += (g_i - result) / t++;
	return result;
}

int main() {
	setlocale(LC_ALL, "Russian");

	array<Particle, N> p;
	creating_field(p);
	recording(p, out);

	//цикл OneChain-ов с записью положений всех частиц 
	size_t iter_num = 0;
	vector<double> g;
	do
	{
		int index = rand() % N;
		OneChain(p, index);
		g = rad_func(p);

		cout << "g = ";
		for (auto i : g)
			cout << to_string(i) << " ";
		cout << endl;
		recording(p, out);
		
	} while (iter_num++ != 10 || avarage(g) > 10.0);

	out.close();
	return 0;
}

/*
1. Создаем поле частиц. Частицы располагаются в правильном (шахматном) порядке, между рядами зазоры в ds.
2. *Определяем ближайших*
3. OneChain:
	for(...){
		выбираем направление вверх/вправо,
		находим ближайшую в этом направлении (перебор ближайших по условиям расстояния центров),
		толкаем частицу и перемещаем её в точку столкновения с ближайшей
		перемещаем флажок
	}
*/
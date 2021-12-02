using namespace std;
#include "Particle.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <vector>
#include <string>
#include <list>
#include <chrono>

#include <fstream>
#include <math.h>
double get_distance_y(const Particle& A, const Particle& B) {
	double dx = B.x - A.x;
	if (dx < 0) dx += L;
	if (dx > L / 2) dx = L - dx;

	double dy = B.y - A.y;
	if (dy < 0) dy += L;

	// check if particle is reachable
	if (dx >= 2 * r) return numeric_limits<double>::max();
	double result = dy - sqrt(4 * r * r - dx * dx);
	return result;
}
double get_distance_x(const Particle& A, const Particle& B) {
	double dy = B.y - A.y;
	if (dy < 0) dy += L;
	if (dy > L / 2) dy = L - dy;

	double dx = B.x - A.x;
	if (dx < 0) dx += L;

	// check if particle is reachable
	if (dy >= 2 * r) return numeric_limits<double>::max();
	double result = dx - sqrt(4 * r * r - dy * dy);
	return result;
}
void creating_field(vector<Particle>& p) {
	int floor = 0, index = 0;
	for (double j = r; j < L; j += (2 * r + ds)) {
		for (double i = r; i < L; i += (2 * r + ds)) {
			p[index].x = i + floor % 2 * r;
			p[index].y = j;
			index++;
		}
		floor++;
	}
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
bool is_reachable(const Particle& a, const Particle& b, bool dir)
{
	if (dir) return get_distance_y(a, b) != numeric_limits<double>::max();
	else return get_distance_x(a, b) != numeric_limits<double>::max();
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
Particle* collide_y(Particle* a, vector<Particle>& A, double& distance)
{
	list<Particle*> reachable_particles;
	for (size_t i = 0; i != N; i++)
		if (&A[i] != a && abs(a->x - A[i].x) < 2 * r)
			reachable_particles.push_back(&A[i]);
	if (reachable_particles.empty()) return a;

	reachable_particles.sort([a](Particle* lhs, Particle* rhs)
		{ return get_distance_y(*a, *lhs) < get_distance_y(*a, *rhs); }
	);

	Particle* b = reachable_particles.front();
	double delta = get_distance_y(*a, *b);
	collision(true, a, delta);
	distance += delta;
	return b;
}
Particle* collide_x(Particle* a, vector<Particle>& A, double& distance)
{
	list<Particle*> reachable_particles;
	for (size_t i = 0; i != N; i++)
		if (&A[i] != a && abs(a->y - A[i].y) < 2 * r)
			reachable_particles.push_back(&A[i]);
	if (reachable_particles.empty()) return a;

	reachable_particles.sort([a](Particle* lhs, Particle* rhs)
		{ return get_distance_x(*a, *lhs) < get_distance_x(*a, *rhs); }
	);

	Particle* b = reachable_particles.front();
	double delta = get_distance_x(*a, *b);
	collision(false, a, delta);
	distance += delta;
	return b;
}
ofstream out("particles.txt", ios::out); // открываем файл для записи
size_t num_recordings = 0;
void recording(vector<Particle>& A, ofstream& out) {
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
	num_recordings++;
}
void OneChain(vector<Particle>& A, int index = 0) { //dir - direction: 0 - right, 1 - up
													   //index - index of start particle, =0
	double distance = 0; //сколько всего пройдено
	bool direction = rand() % 2;
	Particle* a = &A[index];
	while (distance < MAX_ONECHAIN_DISTANCE) {
		double old_distance = distance;
		if (direction) a = collide_y(a, A, distance);
		else a = collide_x(a, A, distance);
		//recording(A, out);
		if (abs(old_distance - distance) < 10e-2) break;
	}

}
vector<double> rad_func(const vector<Particle>& p)
{
	//histo
	vector<int> hist;
	double dr = 0.01 * r;
	for (double d = 0; d <= 0.1; d += dr / r)
	{
		hist.push_back(0);
		for (int i = 0; i < N - 1; i++)
		{
			for (int j = i + 1; j < N; j++)
			{
				/* double dx = p[j].x - p[i].x;
				double dy = p[j].y - p[i].y; */
				if (pow(p[j].x - p[i].x,2) + pow(p[j].y - p[i].y, 2) <= pow((2 * r + d * r) + 0.005, 2))
				{
					if (pow(p[j].x - p[i].x, 2) + (pow(p[j].y - p[i].y, 2) >= pow((2 * r + d * r) - 0.005 , 2)))
						hist[d * 100]++;
				}
			}
		}
	}

	for (int j = 0; j < 11; j++)
	{
		cout << 2 + 0.01 * j << ": " << hist[j] << "\n";
	}
	cout << "\n";
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
	vector<Particle> p(N);
	creating_field(p);
	recording(p, out);

	//цикл OneChain-ов с записью положений всех частиц 
	size_t iter_num = 0;
	vector<double> g;
	auto begin_time = chrono::system_clock::now();
	do
	{
		int index = rand() % N;
		OneChain(p, index);

		/*
		g = rad_func(p);
		cout << "g = ";
		for (auto i : g)
			cout << to_string(i) << " ";
		cout << endl;
		*/
		recording(p, out);
		cout << "Iteration " << iter_num << " ended" << endl;

	} while (iter_num++ != 30 || avarage(g) > 10.0);
	g = rad_func(p);
	cout << "g : \n";
	for (int i = 0; i < 11; i++)
		cout << 2 + 0.01 * i << ":     " << g[i] << " \n";
	cout << endl;
	chrono::duration<double> elapsed = chrono::system_clock::now() - begin_time;
	std::cout << "Main loop time: " << elapsed.count() << std::endl;
	std::cout << "Number of recordings: " << num_recordings << std::endl;
	out.close();
	return 0;
}


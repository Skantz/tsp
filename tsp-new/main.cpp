#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>

using namespace std;

#define TIME_LIMIT 1890 * 1000 * 1000
#define SIZE 1000
#define KOPT_DEPTH 1


void fill_distance_matrix(int n, int (*dm)[SIZE]) {
    double pts[SIZE * 2];
    for (int i = 0; i < 2 * n; i++) {
        cin >> pts[i];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double dx = pts[i * 2] - pts[j * 2];
            double dy = pts[i * 2 + 1] - pts[j * 2 + 1];
            dm[i][j] = (int) round(sqrt(pow(dx, 2) + pow(dy, 2)));
        }
    }
}

vector<int> greedy_tour(int n, const int (*dm)[SIZE]) {
    vector<int> tour(n);

    vector<bool> used(n, false);
    used[0] = true;

    tour.push_back(0);

    for (int i = 1; i < n; i++) {
        int best = -1;

        for (int j = 0; j < n; j++) {
            if (!used[j] && (best == -1 || dm[tour[i - 1]][j] < dm[tour[i - 1]][best])) {
                best = j;
            }
        }

        tour[i] = best;
        used[best] = true;
    }

    return tour;
}

int mod(int a, int b) { return a >= 0 ? a % b : (b - abs(a % b)) % b; }

//TOFIX update maps
int node_shift(vector<int> &tour, int n, const int (*dm)[SIZE]) {
    bool improved = true;

    int capital_delta = 0;

    while (improved) {
        improved = false;
        for (int i = 1; i < n; ++i) {
            for (int j = i + 2; j < n - 2; ++j) {

                //int u = mod((i - 1), n);
                //int v = (j + 1) % n;

                int gain_rand = dm[tour[i - 1]][tour[i]] + dm[tour[i]][tour[i + 1]] + dm[tour[j]][tour[j + 1]];
                int loss_rand = dm[tour[i - 1]][tour[i + 1]] + dm[tour[i]][tour[j]] + dm[tour[i]][tour[j + 1]];

                int delta = (loss_rand) - (gain_rand);

                if (delta < 0) {
                    tour.insert(tour.begin() + j + 1, tour[i]);
                    tour.erase(tour.begin() + i);

                    improved = true;

                    capital_delta += delta;
                }
            }
        }
    }

    return capital_delta;
}

int two_opt(vector<int> &tour, int n, const int (*dm)[SIZE]) {
    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n - 1; j++) {

                int i_m = mod((i - 1), n);
                int j_p = (j + 1) % n;

                int new_cost = dm[tour[i]][tour[j_p]] + dm[tour[i_m]][tour[j]];
                int old_cost = dm[tour[i_m]][tour[i]] + dm[tour[j]][tour[j_p]];

                int delta = new_cost - old_cost;

                if (delta < 0) {
                    reverse(tour.begin() + i, tour.begin() + j + 1);
                    improved = true;

                    capital_delta += delta;
                }
            }
        }
    }

    return capital_delta;
}

int two_opt_move(vector<int> &tour, int n, const int (*dm)[SIZE], int i, int j) {
    int i_m = mod((i - 1), n);
    int j_p = (j + 1) % n;

    int new_cost = dm[tour[i]][tour[j_p]] + dm[tour[i_m]][tour[j]];
    int old_cost = dm[tour[i_m]][tour[i]] + dm[tour[j]][tour[j_p]];

    int delta = new_cost - old_cost;

    reverse(tour.begin() + i, tour.begin() + j + 1);

    return delta;
}

int three_opt_move(vector<int> &tour, int n, const int (*dm)[SIZE], int i, int j, int k) {


    int a, b, c, d, e, f;
    a = tour[i - 1];
    b = tour[i];
    c = tour[j - 1];
    d = tour[j];
    e = tour[k - 1];
    f = tour[k];

    int d0 = dm[a][b] + dm[c][d] + dm[e][f];
    int d1 = dm[a][c] + dm[b][d] + dm[e][f];
    int d2 = dm[a][b] + dm[c][e] + dm[d][f];
    int d3 = dm[a][d] + dm[e][b] + dm[c][f];
    int d4 = dm[f][b] + dm[c][d] + dm[e][a];

    int case_1 = -d0 + d1;
    int case_2 = -d0 + d2;
    int case_3 = -d0 + d3;
    int case_4 = -d0 + d4;
    int case_min = min(case_2, min(case_3, case_4));

    if (case_1 <= case_min) {

        reverse(tour.begin() + i, tour.begin() + j);
        return case_1;

    } else if (case_2 == case_min) {

        reverse(tour.begin() + j, tour.begin() + k);
        return case_2;

    } else if (case_3 == case_min) {

        vector<int> tmp(k - j);

        copy(tour.begin() + j, tour.begin() + k, tmp.begin());
        copy_backward(tour.begin() + i, tour.begin() + j, tour.begin() + k);
        copy(tmp.begin(), tmp.end(), tour.begin() + i);

        return case_3;

    } else {

        reverse(tour.begin() + i, tour.begin() + k);
        return case_4;

    }
}

int tour_distance(vector<int> tour, int n, const int (*dm)[SIZE]) {
    int s = 0;
    for (int i = 0; i < n - 1; i++)
        s += dm[tour[i]][tour[i + 1]];
    return s + dm[tour[n - 1]][tour[0]];
}

vector<int> k_opt(vector<int> tour, int n, const int (*dm)[SIZE]) {
    vector<int> new_tour(tour);

    int delta = 0;
    int u, v, old_dist, new_dist, um, vm, up, vp;
    int a, b, c;
    for (int i = 0; i < KOPT_DEPTH; i++) {

        b = n - 1;
        //TOFIX: allow a == 0

        a = 1 + (int) ((n - 3) * (rand() / (RAND_MAX + 1.0)));
        c = 1 + (int) ((a - 3) * (rand() / (RAND_MAX + 1.0)));
        //a = 0 + (int) ((nbors - 1) * (rand() / (RAND_MAX + 1.0)));
        if (c < 1)
            delta += two_opt_move(new_tour, n, dm, a, b);
        else
            delta += three_opt_move(new_tour, n, dm, c, a, b);

        if (delta < 0)
            break;
    }

    delta += two_opt(new_tour, n, dm);
    delta += node_shift(new_tour, n, dm);

    if (delta < 0)
        return new_tour;

    return tour;
}


int main() {
    srand((unsigned) time(0));

    int n;
    cin >> n;

    int (*dm)[SIZE] = new int[n][SIZE];

    fill_distance_matrix(n, dm);

    vector<int> tour = greedy_tour(n, dm);

    two_opt(tour, n, dm);
    node_shift(tour, n, dm);

    auto start = chrono::high_resolution_clock::now();

    while ((chrono::high_resolution_clock::now() - start).count() < TIME_LIMIT)
        tour = k_opt(tour, n, dm);


    for (auto &itr : tour)
        cout << itr << endl;

    return 0;
}


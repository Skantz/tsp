#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cmath>

using namespace std;

#define TIME_LIMIT 1950 * 1000 * 1000


int tour_distance(vector<int> tour, int n, const vector<vector<int>> &dm) {
    int s = 0;
    for (int i = 0; i < n - 1; i++)
        s += dm[tour[i]][tour[i + 1]];
    return s + dm[tour[n - 1]][tour[0]];
}

inline int mod(int a, int b) { return a >= 0 ? a % b : (b - abs(a % b)) % b; }

inline int randbetween(int a, int b) { return a + (int) ((b - a) * (rand() / (RAND_MAX + 1.0))); }

void brute_force(vector<int> &tour, int n, const vector<vector<int>> &dm) {
    int best = tour_distance(tour, n, dm);

    vector<int> changing_tour(n);

    for (int i = 0; i < n; i++)
        changing_tour[i] = i;

    do {
        int score = tour_distance(changing_tour, n, dm);

        if (score < best) {
            best = score;
            copy(changing_tour.begin(), changing_tour.end(), tour.begin());
        }
    } while (next_permutation(changing_tour.begin() + 1, changing_tour.end()));
}

vector<vector<int>> init_distance_matrix(int n) {
    vector<double> pts(n * 2);
    vector<vector<int>> dm(n, vector<int>(n));

    for (int i = 0; i < n; i++) {
        cin >> pts[2 * i] >> pts[2 * i + 1];

        for (int j = 0; j < i; j++)
            dm[i][j] = dm[j][i] = (int) round(sqrt(
                    pow(pts[i * 2] - pts[j * 2], 2) +
                    pow(pts[i * 2 + 1] - pts[j * 2 + 1], 2)
            ));
    }

    return dm;
}

vector<int> greedy_tour(int n, const vector<vector<int>> &dm) {
    vector<int> tour(n);

    vector<bool> used(n, false);
    used[0] = true;

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

//TOFIX update maps
int node_shift(vector<int> &tour, int n, const vector<vector<int>> &dm) {
    bool improved = true;

    int capital_delta = 0;

    while (improved) {
        improved = false;

        for (int i = 0; i < n - 1; i++) {
            int i_m = mod(i - 1, n);
            int i_p = i + 1;

            for (int j = i + 2; j < n; j++) {
                int j_m = j - 1;

                int delta = dm[tour[i_m]][tour[i_p]] + dm[tour[j_m]][tour[i]] + dm[tour[i]][tour[j]] -
                            (dm[tour[i_m]][tour[i]] + dm[tour[i]][tour[i_p]] + dm[tour[j_m]][tour[j]]);

                if (delta < 0) {
                    tour.insert(tour.begin() + j, tour[i]);
                    tour.erase(tour.begin() + i);

                    improved = true;

                    capital_delta += delta;
                }
            }
        }
    }

    return capital_delta;
}

int two_opt(vector<int> &tour, int n, const vector<vector<int>> &dm) {
    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;

        for (int i = 0; i < n; i++) {
            int i_m = mod(i - 1, n);

            for (int j = i + 2; j < n; j++) {
                int j_m = j - 1;

                int delta = dm[tour[i]][tour[j]] + dm[tour[i_m]][tour[j_m]] -
                            (dm[tour[i_m]][tour[i]] + dm[tour[j_m]][tour[j]]);

                if (delta < 0) {
                    reverse(tour.begin() + i, tour.begin() + j);
                    improved = true;

                    capital_delta += delta;
                }
            }
        }
    }

    return capital_delta;
}

int random_three_opt(vector<int> &tour, int n, const vector<vector<int>> &dm) {
    int i = randbetween(0, n - 4);
    int j = randbetween(i + 2, n - 2);
    int k = randbetween(j + 2, n);

    int a = tour[mod(i - 1, n)];
    int b = tour[i];
    int c = tour[j - 1];
    int d = tour[j];
    int e = tour[k - 1];
    int f = tour[k];

    int d0 = dm[a][b] + dm[c][d] + dm[e][f];
    int d1 = dm[a][c] + dm[b][d] + dm[e][f];
    int d2 = dm[a][b] + dm[c][e] + dm[d][f];
    int d3 = dm[a][d] + dm[e][b] + dm[c][f];
    int d4 = dm[f][b] + dm[c][d] + dm[e][a];
    int d_max = max(d2, max(d3, d4));

    if (d1 <= d_max) {

        reverse(tour.begin() + i, tour.begin() + j);
        return d1 - d0;

    } else if (d2 == d_max) {

        reverse(tour.begin() + j, tour.begin() + k);
        return d2 - d0;

    } else if (d3 == d_max) {

        vector<int> tmp(k - j);

        copy(tour.begin() + j, tour.begin() + k, tmp.begin());
        copy_backward(tour.begin() + i, tour.begin() + j, tour.begin() + k);
        copy(tmp.begin(), tmp.end(), tour.begin() + i);

        return d3 - d0;

    } else {

        reverse(tour.begin() + i, tour.begin() + k);
        return d4 - d0;

    }
}

vector<int> k_opt(vector<int> tour, int n, const vector<vector<int>> &dm) {
    vector<int> new_tour(tour);

    int delta = 0;

    delta += random_three_opt(new_tour, n, dm);
    delta += two_opt(new_tour, n, dm);
    delta += node_shift(new_tour, n, dm);

    if (delta < 0)
        return new_tour;

    return tour;
}


int main() {
    auto start = chrono::high_resolution_clock::now();

    srand((unsigned) time(nullptr));

    int n;
    cin >> n;

    vector<vector<int>> dm = init_distance_matrix(n);

    vector<int> tour = greedy_tour(n, dm);

    if (n < 10) {
        brute_force(tour, n, dm);
    } else {

        two_opt(tour, n, dm);
        node_shift(tour, n, dm);

        while ((chrono::high_resolution_clock::now() - start).count() < TIME_LIMIT)
            tour = k_opt(tour, n, dm);
    }

    for (auto &itr : tour)
        cout << itr << endl;

    return 0;
}


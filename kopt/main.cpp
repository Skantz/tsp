#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include <algorithm>
#include <map>
#include <random>
#include <ctime>
#include <chrono>

using namespace std;

#define TIME_LIMIT 1890 * 1000 * 1000

int n;
vector<vector<int>> dm;
vector<int> heuristic_tour;
int heuristic_cost;
vector<vector<int>> neighbours;
chrono::time_point<chrono::high_resolution_clock> start_time;

int mod(int a, int b) { return a >= 0 ? a % b : (b - abs(a % b)) % b; }

class Tour {
    vector<int> tour;
    set<pair<int, int>> edges;
    int size;

public:
    explicit Tour(vector<int> &tour) : tour(tour), size(tour.size()) {
        int u = tour[0], v = tour[size - 1];

        edges.insert(u < v ? pair(u, v) : pair(v, u));

        for (int i = 1; i < size; i++) {
            u = tour[i], v = tour[i - 1];
            edges.insert(u < v ? pair(u, v) : pair(v, u));
        }
    }

    int at(int i) {
        return tour[i];
    }

    bool contains(pair<int, int> p) {
        return edges.find(p) != edges.end();
    }

    bool contains(int u, int v) {
        return edges.find(u < v ? pair(u, v) : pair(v, u)) != edges.end();
    }

    int index(int node) {
        for (int i = 0; i < size; i++)
            if (tour[i] == node)
                return i;

        return -1;
    }

    vector<int> around(int node) {
        for (int i = 0; i < size; i++)
            if (tour[i] == node)
                return vector<int>{tour[mod(i - 1, size)], tour[(i + 1) % size]};

        return vector<int>{};
    }

    int pred(int i) {
        return tour[mod(i - 1, size)];
    }

    int succ(int i) {
        return tour[(i + 1) % size];
    }

    pair<bool, vector<int>> generate(const set<pair<int, int>> &broken, const set<pair<int, int>> &joined) {
        set<pair<int, int>> tmp;
        set<pair<int, int>> new_edges;
        set_difference(edges.begin(), edges.end(), broken.begin(), broken.end(), inserter(tmp, tmp.end()));
        set_union(tmp.begin(), tmp.end(), joined.begin(), joined.end(), inserter(new_edges, new_edges.end()));

        if (new_edges.size() < size)
            return pair(false, vector<int>());

        map<int, int> successors;
        int node = 0;

//        for (auto & edge : new_edges) {
//            cout << edge.first << " " << edge.second << " - ";
//        }
//        cout << endl;

        towhile:
        if (!new_edges.empty()) {
            for (auto &edge : new_edges) {
                if (edge.first == node) {
                    successors.insert_or_assign(node, edge.second);
                    node = edge.second;
                    new_edges.erase(edge);
                    goto towhile;
                }
                if (edge.second == node) {
                    successors.insert_or_assign(node, edge.first);
                    node = edge.first;
                    new_edges.erase(edge);
                    goto towhile;
                }
            }

            new_edges.clear();
        }

//        cout << new_edges.size() << " " << successors.size() << endl;

        if (successors.size() < size)
            return pair(false, vector<int>());

        int succs = successors[0];
        vector<int> new_tour;
        new_tour.push_back(0);
        set<int> visited = {0};

        while (visited.find(succs) == visited.end()) {
            visited.insert(succs);
            new_tour.push_back(succs);
            succs = successors[succs];
        }

        return pair(new_tour.size() == size, new_tour);
    }
};

int tour_cost(const vector<int> &tour) {
    int cost = dm[tour[tour.size() - 1]][tour[0]];

    for (int i = 1; i < tour.size(); i++)
        cost += dm[tour[i - 1]][tour[i]];

    return cost;
}

void random_distance_matrix() {
    srand((unsigned) time(0));

    n = 10;
    int spread = 100;

    vector<double> pts(n * 2);
    for (int i = 0; i < 2 * n; i++) {
        pts[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/spread));
    }

    dm = vector(n, vector<int>(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double dx = pts[i * 2] - pts[j * 2];
            double dy = pts[i * 2 + 1] - pts[j * 2 + 1];
            dm[i][j] = (int) round(sqrt(pow(dx, 2) + pow(dy, 2)));
        }
    }
}

void init_distance_matrix() {
    cin >> n;

    vector<double> pts(n * 2);
    for (int i = 0; i < 2 * n; i++) {
        cin >> pts[i];
    }

    dm = vector(n, vector<int>(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double dx = pts[i * 2] - pts[j * 2];
            double dy = pts[i * 2 + 1] - pts[j * 2 + 1];
            dm[i][j] = (int) round(sqrt(pow(dx, 2) + pow(dy, 2)));
        }
    }
}

vector<int> greedy_tour() {
    vector<int> tour(n);

    vector<bool> used(n, false);
    used[0] = true;
    tour[0] = 0;

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


vector<pair<int, pair<int, int>>>
closest(int t2i, Tour &tour, int gain, const set<pair<int, int>> &broken, const set<pair<int, int>> &joined) {
    map<int, pair<int, int>> nbors;

    for (auto &node : neighbours[t2i]) {
        auto yi = t2i < node ? pair(t2i, node) : pair(node, t2i);
        int Gi = gain - dm[t2i][node];

        if (Gi <= 0 || broken.find(yi) != broken.end() || tour.contains(yi))
            continue;

        for (auto &succ : tour.around(node)) {
            pair<int, int> xi = node < succ ? pair(node, succ) : pair(succ, node);

            if (broken.find(xi) == broken.end() && joined.find(xi) == joined.end()) {
                int diff = dm[node][succ] - dm[t2i][node];

                if (nbors.find(node) != nbors.end() && diff > nbors[node].first)
                    nbors[node].first = diff;
                else
                    nbors[node] = pair(diff, Gi);
            }
        }
    }

    // FIXME: I hate this
    vector<pair<int, pair<int, int>>> nborsvec(nbors.size());

    copy(nbors.begin(), nbors.end(), nborsvec.begin());

    sort(nborsvec.begin(), nborsvec.end(),
         [](const pair<int, pair<int, int>> &x, const pair<int, pair<int, int>> &y) {
             return x.second.first > y.second.first;
         });

    return nborsvec;
}


bool
choose_y(Tour &tour, int t1, int t2i, int gain, const set<pair<int, int>> &broken, const set<pair<int, int>> &joined);

bool
choose_x(Tour &tour, int t1, int last, int gain, const set<pair<int, int>> &broken, const set<pair<int, int>> &joined) {
    auto around = tour.around(last);

    if (broken.size() == 4) {
        if (dm[around[0]][last] > dm[around[1]][last])
            around.pop_back();
        else
            around.erase(around.begin());
    }

    for (int &t2i : around) {
        if ((chrono::high_resolution_clock::now() - start_time).count() > TIME_LIMIT)
            return false;

        pair<int, int> xi = last < t2i ? pair(last, t2i) : pair(t2i, last);
        int Gi = gain + dm[last][t2i];

        if (joined.find(xi) == joined.end() && broken.find(xi) == broken.end()) {
            set<pair<int, int>> added;
            copy(joined.begin(), joined.end(), inserter(added, added.end()));
            set<pair<int, int>> removed;
            copy(broken.begin(), broken.end(), inserter(removed, removed.end()));


            removed.insert(xi);
            added.insert(t2i < t1 ? pair(t2i, t1) : pair(t1, t2i));

            int relink = Gi - dm[t2i][t1];
            pair<bool, vector<int>> is_new_tour = tour.generate(removed, added);

            if (!is_new_tour.first && added.size() > 2)
                continue;

//            printf("Choose x: (%d, %d) = %d\n", xi.first, xi.second, dm[last][t2i]);

            if (is_new_tour.first && relink > 0) {
                heuristic_tour = is_new_tour.second;
                heuristic_cost -= relink;

//                printf("Found %d-opt with gain %d\n", added.size(), relink);

                return true;
            } else {
                bool choice = choose_y(tour, t1, t2i, Gi, removed, joined);

                if (broken.size() == 2 && choice)
                    return true;

                return choice;
            }
        }
    }

    return false;
}

bool
choose_y(Tour &tour, int t1, int t2i, int gain, const set<pair<int, int>> &broken, const set<pair<int, int>> &joined) {
    auto close = closest(t2i, tour, gain, broken, joined);

    int top = broken.size() == 2 ? 5 : 1;

    for (pair<int, pair<int, int>> &it : close) {
        if ((chrono::high_resolution_clock::now() - start_time).count() > TIME_LIMIT)
            return false;

        int node = it.first;
        int Gi = it.second.second;

        pair<int, int> yi = t2i < node ? pair(t2i, node) : pair(node, t2i);

        set<pair<int, int>> added;
        copy(joined.begin(), joined.end(), inserter(added, added.end()));
        added.insert(yi);

//        printf("Choose y: (%d, %d) = %d\n", yi.first, yi.second, dm[t2i][node]);

        if (choose_x(tour, t1, node, Gi, broken, added))
            return true;

        if (--top == 0)
            return false;
    }

    return false;
}


bool improve() {
    Tour tour = Tour(heuristic_tour);

    for (auto &t1 : heuristic_tour) {
        auto around = tour.around(t1);

        for (auto &t2: around) {
            set<pair<int, int>> broken = {t1 < t2 ? pair(t1, t2) : pair(t2, t1)};
            int gain = dm[t1][t2];
            auto close = closest(t2, tour, gain, broken, set<pair<int, int>>());

            int tries = 5;

            for (auto &it: close) {
                if ((chrono::high_resolution_clock::now() - start_time).count() > TIME_LIMIT)
                    return false;

                int t3 = it.first;
                int Gi = it.second.second;

                if (t3 == around[0] || t3 == around[1])
                    continue;

                set<pair<int, int>> joined = {t2 < t3 ? pair(t2, t3) : pair(t3, t2)};

//                printf("Start: X: (%d, %d) - Y: (%d, %d) = %d\n", t1, t2, t2, t3, Gi);

                if (choose_x(tour, t1, t3, Gi, broken, joined))
                    return true;

                if (--tries == 0)
                    break;
            }
        }
    }

    return false;

}

void optimise() {
    bool better = true;

    // FIXME: do something with the 'self.solutions'

    neighbours = vector<vector<int>>(n, vector<int>());

    for (int &i : heuristic_tour)
        for (int &j : heuristic_tour)
            if (dm[i][j] > 0)
                neighbours[i].push_back(j);

    while (better && (chrono::high_resolution_clock::now() - start_time).count() < TIME_LIMIT) {
        better = improve();
        // FIXME: 'self.solutions'
    }
}

void print_result() {
    for (int &node : heuristic_tour)
        cout << node << endl;
}


int main() {
    start_time = chrono::high_resolution_clock::now();

//    random_distance_matrix();
    init_distance_matrix();

    heuristic_tour = greedy_tour();
    heuristic_cost = tour_cost(heuristic_tour);

//    for (int &node : heuristic_tour)
//        cout << node << " ";
//
//    cout << endl << heuristic_cost << endl;

    optimise();

//    for (int &node : heuristic_tour)
//        cout << node << " ";
//
//    cout << endl << heuristic_cost << endl;

    print_result();

    return 0;
}

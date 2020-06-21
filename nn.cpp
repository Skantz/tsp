#include <iostream>
#include <sstream>
#include <istream>
#include <vector>
#include <math.h> 
#include <algorithm>
#include <assert.h>
#include <ctime>
#include <chrono>
#include <functional>
#include <queue>
#include <map>


using namespace std;

#define DEBUG 0 
#define VERBOSE 1



//x1, y1, x2, y2..
vector<vector<int>> create_distance_matrix(vector<double> points) {
    int n = points.size() / 2;
    vector<vector<int>> dm;
    for (int i = 0; i < n ; i++) {
        dm.push_back(vector<int>());
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) { 
            double x_dist = (points[i*2] - points[j*2]);
            double y_dist = (points[i*2 + 1] - points[j*2 + 1]);
            dm[i].push_back(round(sqrt((pow(x_dist, 2) + pow(y_dist, 2)))));
        }
    }
    return dm;
}

vector<vector<int>> create_nearest_neighbor_matrix(vector<vector<int>> dm, int k) {
    int n = dm[0].size();
    vector<vector<int>> nm;
    for (int i = 0; i < n ; i++) {
        nm.push_back(vector<int>(0, k));
    }

    for (int r = 0; r < n; r++) {

        vector<int> row = dm[r];
        priority_queue<pair<int, int>> row_queue;

        for (int i = 0; i < row.size(); ++i) {
            row_queue.push(std::pair<double, int>(row[i], i));
        }

        for (int i = 0; i < k; ++i) {
            int k_idx = row_queue.top().second;
            nm[r].push_back(k_idx);
            row_queue.pop();
        }
    }
    return nm;
}

int tour_distance(vector<int> tour, vector<vector<int>> dm) {
    
    int s = 0;
    for (int i = 0; i < tour.size() - 1; i++) {
        s += dm[tour[i]][tour[i + 1]];
    }
    return s + dm[tour[tour.size() - 1]][tour[0]];
}

vector<int> greedy_tour(vector<vector<int>>& dm) {

    int n = dm[0].size();
    vector<int> tour;
    vector<int> used_vertices(n,0);

    int min_pos = 1;
    used_vertices[0] = 1;
    tour.push_back(0);
    for (int i = 1; i < n; i++) {
        min_pos = -1;
        for (int j = 0; j < n; j++) {
        int previous = tour[-1];
        if (used_vertices[j] != 1 && (min_pos == -1 || (dm[previous][j] < dm[previous][min_pos]))) 
            min_pos = j;
        }
    tour.push_back(min_pos);
    used_vertices[min_pos] = 1;
    }

    return tour;
}

int mod (int a, int b)
{ return a >= 0 ? a % b : ( b - abs ( a % b ) ) % b; }

bool is_default_argument(vector<vector<int>>& matrix) {
    if (matrix.size() == 1 && matrix[0].size() == 0 && matrix[0][0] == 0)
        return true;
    return false;
}


int two_opt(vector<int>& tour, vector<vector<int>>& dm, const vector<vector<int>>& nm = vector<vector<int>>{{0}}) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;
        for (int i = 0; i < tour.size() - 1; i++) {
            for (int j = i + 1; j < tour.size() - 1; j++) {

                int i_m = mod((i - 1), n);
                int i_p = mod((i + 1), n);
                int j_m = mod((j - 1), n);
                int j_p = mod((j + 1), n);  
                //cout << "i_ m " << i_m << " j p" << j_p << endl;
                int cost_1 = dm[tour[i_m]][tour[i]];
                int cost_2 = dm[tour[j]][tour[j_p]];

                int new_cost_1 = dm[tour[i]][tour[j_p]];
                int new_cost_2 = dm[tour[i_m]][tour[j]];

                int delta = (new_cost_1 + new_cost_2) - (cost_1 + cost_2);

                if (delta < 0) {
                    //SWAP
                    #if DEBUG == 1
                    double old_cost = tour_distance(tour, dm);
                    #endif
                    std::reverse(tour.begin() + i , tour.begin() + j + 1);
                    improved = true;
                    #if DEBUG == 1
                    double new_cost = tour_distance(tour, dm);
                    //cout << old_cost << " -> " << new_cost << endl;
                    assert(new_cost < old_cost);
                    assert(new_cost == old_cost + delta);
                    #endif 
                    capital_delta += delta;
                }
            }
        }
    }
    return capital_delta;
}


int two_opt_nn(vector<int>& tour, vector<vector<int>>& dm, const vector<vector<int>>& nm, map<int, int>& idx_map) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;
        for (int i = 0; i < tour.size() - 1; i++) {
            int p1 = tour[i];
            for (auto& p2: nm[p1]) {
                for (int b = 0; b < 2; b++) {
                    int j = idx_map[p2] - 1;
                    if (b != 0) {
                        int j = idx_map[p2] + 1;
                    }

                    //TOFIX: support j == n - 1
                    if (i >= j || j == n - 1) {
                        continue;
                    }
                    //cout << "i " << tour[idx_map[p1]] << " -> " << idx_map[p1] << " -> " <<  p1 << " -> " << i << endl;
                    //cout << "j " << tour[idx_map[p2]] << " -> " << idx_map[p2] << " -> " <<  p2 << " -> " << j << endl;
                    assert(tour[idx_map[p1]] == p1);
                    assert(tour[idx_map[p2]] == p2);
                    int i_m = mod((i - 1), n);
                    int i_p = mod((i + 1), n);
                    int j_m = mod((j - 1), n);
                    int j_p = mod((j + 1), n);  
                    //cout << "i_ m " << i_m << " j p" << j_p << endl;
                    int cost_1 = dm[tour[i_m]][tour[i]];
                    int cost_2 = dm[tour[j]][tour[j_p]];

                    int new_cost_1 = dm[tour[i]][tour[j_p]];
                    int new_cost_2 = dm[tour[i_m]][tour[j]];

                    int delta = (new_cost_1 + new_cost_2) - (cost_1 + cost_2);

                    if (delta < 0) {
                        //SWAP
                        #if DEBUG == 1
                        double old_cost = tour_distance(tour, dm);
                        #endif
                        //[i, j+1)
                        assert(i <= j);
                        for (int t = 0; t < (j - i) + 1; t++ ) {
                            if (i + t > j - t) break;
                            int s = idx_map[tour[i + t]];
                            idx_map[tour[i + t]] = idx_map[tour[j - t]];
                            idx_map[tour[j - t]] = s;
                        }
                        std::reverse(tour.begin() + i , tour.begin() + j + 1);
                        improved = true;
                        #if DEBUG == 1
                        double new_cost = tour_distance(tour, dm);
                        //cout << old_cost << " -> " << new_cost << endl;
                        assert(new_cost < old_cost);
                        assert(new_cost == old_cost + delta);
                        #endif 
                        capital_delta += delta;
                        //cout << tour[idx_map[p1]] << " -> " << idx_map[p1] << " -> " <<  p1 << " -> " << i << endl;
                        // p1 = tour[i];

                        int si = i;
                        //idx_map[p1] = j;
                        #if DEBUG == 1
                        cout << "i -> j " << i << " -> " << j << endl; 
                        #endif
                        //idx_map[p2] = i ;
                        //assert(tour[j] == p1);
                        //assert(tour[i] == p2);
                        break;
                    }
                }
            }
            if (improved) {
                break;
            }
        }
    }
    return capital_delta;
}


int two_opt_move(vector<int>& tour, vector<vector<int>> dm, map<int, int>& idx_map, int i, int j) {

    #if DEBUG == 1
    vector<int> saved_tour(tour);
    int p1 = tour[i];
    int p2 = tour[j];
    #endif

    int n = tour.size();
    
    int i_m = mod((i - 1), n);
    int i_p = mod((i + 1), n);
    int j_m = mod((j - 1), n);
    int j_p = mod((j + 1), n);  

    int cost_1 = dm[tour[i_m]][tour[i]];
    int cost_2 = dm[tour[j]][tour[j_p]];

    int new_cost_1 = dm[tour[i]][tour[j_p]];
    int new_cost_2 = dm[tour[i_m]][tour[j]];

    int delta = (new_cost_1 + new_cost_2) - (cost_1 + cost_2);

    for (int t = 0; t < (j - i) + 1; t++ ) {
        if (i + t > j - t) break;
        int s = idx_map[tour[i + t]];
        idx_map[tour[i + t]] = idx_map[tour[j - t]];
        idx_map[tour[j - t]] = s;
    }
    std::reverse(tour.begin() + i , tour.begin() + j + 1);

    #if DEBUG == 1
    assert(tour_distance(tour, dm) ==  delta + tour_distance(saved_tour, dm));
    if (delta < 0) assert(tour_distance(tour, dm) < tour_distance(saved_tour, dm));
    assert(tour[j] == p1);
    assert(tour[i] == p2);
    assert(tour[idx_map[p1]] == p1);
    assert(tour[idx_map[p2]] == p2);
    #endif

    return delta;
}

vector<int> k_opt(vector<int> tour, vector<vector<int>> dm, vector<vector<int>> nm, map<int, int>& idx_map, int depth) {
    // 4 3 2 scheme

    int ts = tour.size();
    assert(ts >= 1);
    vector<int> new_tour(tour);
    map<int, int> saved_map;
        for (int i = 0; i < ts; i++) {
            saved_map[i] = idx_map[i];
        }
    int delta = 0;
    int a = (int) (ts * (rand() / (RAND_MAX + 1.0)));
    //int b = (int) (ts * (rand() / (RAND_MAX + 1.0)));
    int b = ts - 1;

    for (int i = 0; i < depth; i++) {
        //a = b;
        b = ts - 1;
        //TOFIX: allow a == 0
        a = 1 + (int) ((ts - 1) * (rand() / (RAND_MAX + 1.0)));
        //b = mod(b + 1, n);
        a = min(a, b);
        b = max(a, b);   

        delta += two_opt_move(new_tour, dm, idx_map, a, b);
    }

    delta += two_opt_nn(new_tour, dm, nm, idx_map);
    //cout << tour_distance(tour, dm) << " " <<  tour_distance(saved_tour, dm) << " " << delta;
    //cout << "Kopt delta " << delta << endl;
    #if DEBUG == 1
    vector<int> saved_tour(tour);
    #endif
    if (delta < 0) {
        tour = new_tour;
        #if DEBUG == 1
        cout << "delta " << delta << endl;
        assert(tour_distance(tour, dm) < tour_distance(saved_tour, dm));
        #endif
    }
    else if (delta > 0) {
        for (int i = 0; i < ts; i++) {
            idx_map[i] = saved_map[i];
        }
        #if DEBUG == 1
        cout << "delta " << delta << endl;
        assert(tour_distance(new_tour, dm) > tour_distance(saved_tour, dm));
        #endif
    }
    return tour;

}


int main(int argc, char *argv[]) {

    srand( (unsigned)time(0) );
    vector<double> points;

    int n;
    cin >> n;

    double val;
    for (int i = 0; i < 2*n; i++) {
        std::cin >> val;
        points.push_back(val);
    }

    vector<vector<int>> dm = create_distance_matrix(points);
    vector<vector<int>> nm = create_nearest_neighbor_matrix(dm, 5);
    
    vector<int> tour;
    map<int, int> place_to_index_map;
    for (int i = 0; i < n; i++) {
        tour.push_back(i);
        place_to_index_map[i] = i;
    }

    //random_shuffle (tour.begin(), tour.end());
    //if (tour[0] % 2 == 0) {
    //    tour = greedy_tour(dm);
    //}

    two_opt(tour, dm);
    //two_opt_nn(tour, dm, nm, place_to_index_map);
    for (int j = 0; j < n; j++) {
        place_to_index_map[j] = distance(tour.begin(),  find(tour.begin(), tour.end(), j));
    }



    auto t1 = chrono::high_resolution_clock::now();
    auto t2 = chrono::high_resolution_clock::now();

    #if VERBOSE == 1
    cout << "Finish 2-opt" << endl;
    int best_tour_dist = *max_element(dm[0].begin(), dm[0].end()) * n  + 1;
    #endif

    int non_improvement_threshold = 10;
    long raise_depth_every = 1000;
    int depth = 3;
    int i_;

    for (int i = 0; i < 5000000; i++) {
        //int depth = n <= 200? 10: n <= 500? 4: n <= 750? 3 : 2;
        tour = k_opt(tour, dm, nm, place_to_index_map, depth);
        
        //TOFIX: map should not be invalid between calls to k_opt
        for (int j = 0; j < n; j++) {
            place_to_index_map[j] = distance(tour.begin(),  find(tour.begin(), tour.end(), j));
        }

        t2 = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> time_since = t2 - t1;
        #if VERBOSE == 1
        int new_dist = tour_distance(tour, dm);
        if (new_dist < best_tour_dist) {
            cout << new_dist << " <- " << best_tour_dist << " in " << time_since.count() << endl;
            best_tour_dist = new_dist;
        }
        #endif
        #if DEBUG == 1
        for (int i = 0; i < n; i++) {
            cout << i << " " << tour[place_to_index_map[i]] << endl;
            assert(i == tour[place_to_index_map[i]]);
        }
        #endif

        if (time_since.count() > 1910) {
            i_ = i;
            break;
        }


    }

    #if VERBOSE == 1
    cout << tour_distance(tour, dm) << endl;
    cout << "n steps " << i_ << endl;
    #endif

    for (int i = 0; i < n; i++) {
        cout << tour[i];
        if (i < n-1)
            cout << "\n";
    }

  return 0;

}

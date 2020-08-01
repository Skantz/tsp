
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
#include <set>

using namespace std;

#define DEBUG   0
#define VERBOSE 0

struct tour_o {

    vector<int> id;
    vector<int> pos;

    void set(int a, int b) {
        pos[b]  = a;
        id[b]   = a;
    }

    void swap(int a, int b) {
        int was_at_a  = id[a];
        int was_at_b  = id[b];
        pos[was_at_a] = b;
        id[b] = was_at_a;
        pos[was_at_b] = a;
        id[a] = was_at_b;
    }

    void reverse(int a, int b) {
        assert(b >= a);
        
        #if DEBUG == 1
        bool success = false;
        #endif
        //cout  << " a " << a << " b " << b << endl;
        for (int i = 0; i < (b - a)/2 + 1; i++) {
            //cout << "i " << i << " a " << a << " b " << b << endl;
            swap(a + i, b - i);
            #if DEBUG == 1
            if (a + i >= b - i - 1) {success = true;}
            if (a + i > b - i) {success = false; break;}
            #endif

        }   

        #if DEBUG == 1
        assert(success == true);
        #endif
    }

    int size() {
        return id.size();
    }

    int at(int idx) {
        #if DEBUG == 1
        assert(0 <= idx && idx < id.size());
        #endif
        return id[idx];
    }

    tour_o() {
        id = vector<int>();
        pos = vector<int>();
    };

    tour_o(int n) {
        id = vector<int>(n);
        pos = vector<int>(n);
        for (int i = 0; i < n; i++) {
            id[i] = i;
            pos[i] = i;
    }  };

    tour_o(vector<int> v) {
        int n = v.size();
        id = vector<int>(n);
        pos = vector<int>(n);
        for (int i = 0; i < v.size(); i++) {
            id[i] = i;
            pos[i] = v[i];
    }   }

    /*
    tour_o(tour_o& t) {
        for(int i = 0; i < t.id.size(); i++) {
            id[i] = t.id[i];
            pos[i] = t.pos[i];
        }
    };
    */
};

void nothing() {
    ;
}

bool all_different(vector<int> t) {
    vector<int> mem(t.size(), 0);
    for (auto n: t) {
        mem[n] += 1;
        if (mem[n] > 1)
            return false;
    }
    return true;
}

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
    int n = tour.size();
    for (int i = 0; i < n - 1; i++) {
        s += dm[tour[i]][tour[i + 1]];
    }
    return s + dm[tour[n - 1]][tour[0]];
}



int tour_distance(tour_o tour , vector<vector<int>> dm) {
    
    int s = 0;
    int n = tour.size();
    for (int i = 0; i < n - 1; i++) {
        s += dm[tour.at(i)][tour.at(i + 1)];
    }
    return s + dm[tour.at(n - 1)][tour.at(0)];
}



vector<int> greedy_tour(const vector<vector<int>>& dm) {

    int n = dm[0].size();
    vector<int> tour;

    vector<int> used_vertices(n,0);
    
    used_vertices[0] = 1;
    tour.push_back(0);

    for (int i = 1; i < n; i++) {
        int min_pos = -1;
        for (int j = 0; j < n; j++) {
            int previous = tour[i-1];
            if (used_vertices[j] != 1 && (min_pos == -1 || (dm[previous][j] < dm[previous][min_pos]))) {
                min_pos = j;
            }
        }
    tour.push_back(min_pos);
    used_vertices[min_pos] = 1;
    }

    return tour;
}


void greedy_tour(tour_o* tour, const vector<vector<int>>& dm) {

    int n = dm[0].size();

    vector<int> used_vertices(n,0);

    used_vertices[0] = 1;
    (*tour).set(0, 0);

    for (int i = 1; i < n; i++) {
        int min_pos = -1;
        for (int j = 0; j < n; j++) {
            int previous = (*tour).at(i-1);
            if (used_vertices[j] != 1 && (min_pos == -1 || (dm[previous][j] < dm[previous][min_pos]))) {
                min_pos = j;
            }
        }
    (*tour).set(i, min_pos);
    used_vertices[min_pos] = 1;
    }

    return;
}


int mod (int a, int b)
{ return a >= 0 ? a % b : ( b - abs ( a % b ) ) % b; }


//TOFIX update maps
int node_shift(vector<int>& tour, const vector<vector<int>>& dm) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;
        for (int i = 1; i < n ; ++i) {
            for (int j = i + 2; j < n - 2; ++j) {

                //int v1 = mod((i - 1), n);
                //int u1 = (j + 1) % n; 
                int v1 = i;
                int u1 = j;

                int gain_rand = dm[tour[v1 - 1]][tour[v1]]  +  dm[tour[v1]][tour[v1 + 1]]  +  dm[tour[u1]][tour[u1+1]];
                int loss_rand = dm[tour[v1 - 1]][tour[v1 + 1]]  +  dm[tour[v1]][tour[u1]]  +  dm[tour[v1]][tour[u1 + 1]];

                int delta = (loss_rand) - (gain_rand);

                if (delta < 0) {

                    #if DEBUG == 1
                    double old_cost = tour_distance(tour, dm);
                    #endif
                    /*
                    int saved_i = tour[i];
                    tour[i] = tour[j];
                    tour[j] = saved_i;
                    */
                    tour.insert(tour.begin() + u1 + 1, tour[i]);
                    tour.erase(tour.begin() + v1);

                    improved = true;

                    #if DEBUG == 1
                    double new_cost = tour_distance(tour, dm);
                    //cout << old_cost << " -> " << new_cost <<  " (node shift)" << endl;
                    assert(new_cost == old_cost + delta);
                    assert(new_cost < old_cost);
                    #endif 
                    
                    capital_delta += delta;
    }   }   }   }

    return capital_delta;
}

int two_opt(vector<int>& tour, const vector<vector<int>>& dm) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;
        for (int i = 0; i < n ; ++i) {
            for (int j = i + 1; j < n - 1; ++j) {

                int i_m = mod((i - 1), n);
                int j_p = (j + 1) % n; 

                int cost_1 = dm[tour[i_m]][tour[i]];
                int cost_2 = dm[tour[j]][tour[j_p]];

                int new_cost_1 = dm[tour[i]][tour[j_p]];
                int new_cost_2 = dm[tour[i_m]][tour[j]];

                int delta = (new_cost_1 + new_cost_2) - (cost_1 + cost_2);

                if (delta < 0) {

                    #if DEBUG == 1
                    double old_cost = tour_distance(tour, dm);
                    #endif

                    std::reverse(tour.begin() + i, tour.begin() + j + 1);
                    improved = true;

                    #if DEBUG == 1
                    double new_cost = tour_distance(tour, dm);
                    //cout << old_cost << " -> " << new_cost << endl;
                    assert(new_cost < old_cost);
                    assert(new_cost == old_cost + delta);
                    #endif 
                    
                    capital_delta += delta;
    }   }   }   }

    return capital_delta;
}


int two_opt(tour_o& tour, const vector<vector<int>>& dm) {

    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;
        for (int i = 0; i < n ; ++i) {
            for (int j = i + 1; j < n - 1; ++j) {

                int i_m = mod((i - 1), n);
                int j_p = (j + 1) % n; 

                int cost_1 = dm[tour.at(i_m)][tour.at(i)];
                int cost_2 = dm[tour.at(j)][tour.at(j_p)];

                int new_cost_1 = dm[tour.at(i)][tour.at(j_p)];
                int new_cost_2 = dm[tour.at(i_m)][tour.at(j)];

                int delta = (new_cost_1 + new_cost_2) - (cost_1 + cost_2);

                if (delta < 0) {

                    #if DEBUG == 1
                    double old_cost = tour_distance(tour, dm);
                    #endif

                    tour.reverse(i, j);
                    improved = true;

                    #if DEBUG == 1
                    double new_cost = tour_distance(tour, dm);
                    //cout << old_cost << " -> " << new_cost << endl;
                    assert(new_cost < old_cost);
                    assert(new_cost == old_cost + delta);
                    #endif 
                    
                    capital_delta += delta;
    }   }   }   }

    return capital_delta;
}


int two_opt_nn(vector<int>& tour, const vector<vector<int>>& dm, vector<vector<int>> nm, vector<set<int>>& ns) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;
        for (int i = 0; i < n ; ++i) {

            int i_m = mod((i - 1), n);
            int i_p = (i + 1) % n;

            for (auto& i_nbor: nm[i]) {

                for (int direction = -1; direction <= 1; direction += 2) {

                    int j = i_nbor + direction;

                    int j_p = (j + 1) % n; 

                    int cost_1 = dm[tour[i_m]][tour[i]];
                    int cost_2 = dm[tour[j]][tour[j_p]];

                    int new_cost_1 = dm[tour[i]][tour[j_p]];
                    int new_cost_2 = dm[tour[i_m]][tour[j]];

                    int delta = (new_cost_1 + new_cost_2) - (cost_1 + cost_2);

                    if (delta < 0) {

                        #if DEBUG == 1
                        double old_cost = tour_distance(tour, dm);
                        #endif

                        std::reverse(tour.begin() + i, tour.begin() + j + 1);
                        improved = true;

                        #if DEBUG == 1
                        double new_cost = tour_distance(tour, dm);
                        //cout << old_cost << " -> " << new_cost << endl;
                        assert(new_cost < old_cost);
                        assert(new_cost == old_cost + delta);
                        #endif 
                        
                        capital_delta += delta;
    }   }   }   }   }

    return capital_delta;
}



int two_opt_nn(tour_o tour, const vector<vector<int>>& dm, vector<vector<int>> nm, vector<set<int>>& ns) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;
        for (int i = 0; i < n ; ++i) {

            int i_m = mod((i - 1), n);
            int i_p = (i + 1) % n;

            for (auto& i_nbor: nm[i]) {

                for (int direction = -1; direction <= 1; direction += 2) {

                    int j = tour.at(i_nbor) + direction;
                    j = j % n;

                    int j_p = (j + 1) % n; 

                    int cost_1 = dm[tour.at(i_m)][tour.at(i)];
                    int cost_2 = dm[tour.at(j)][tour.at(j_p)];

                    int new_cost_1 = dm[tour.at(i)][tour.at(j_p)];
                    int new_cost_2 = dm[tour.at(i_m)][tour.at(j)];

                    int delta = (new_cost_1 + new_cost_2) - (cost_1 + cost_2);

                    if (delta < 0) {

                        #if DEBUG == 1
                        double old_cost = tour_distance(tour, dm);
                        #endif

                        //std::reverse(tour.begin() + i, tour.begin() + j + 1);
                        tour.reverse(i, j);
                        improved = true;

                        #if DEBUG == 1
                        double new_cost = tour_distance(tour, dm);
                        //cout << old_cost << " -> " << new_cost << endl;
                        assert(new_cost < old_cost);
                        assert(new_cost == old_cost + delta);
                        #endif 
                        
                        capital_delta += delta;
    }   }   }   }   }

    return capital_delta;
}


//TOFIX - Not working! It needs a map from ids to pos
int two_opt_no_look(vector<int>& tour, const vector<vector<int>>& dm) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;

    vector<bool> bits(n, false);
    while (improved) {
        improved = false;
        for (int i = 0; i < n ; i++) {
            if (bits[i] == true) {continue;}
            for (int j = i + 1; j < n - 1; j++) {
                if (bits[j] == true) {continue;}
                int i_m = mod((i - 1), n);
                int j_p = mod((j + 1), n);  

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
                    //i -= 1;
                    //break;
                }
            }
            if (!improved) {
                bits[i] = true;
            }
        }
    }
    return capital_delta;
}


//Not working, indices have many errors.
int or_opt(vector<int>& tour, const vector<vector<int>>& dm) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;
    while (improved) {
        improved = false;
        for (int sl= 3; sl < 4; sl++) {
            //TOFIX: allow i in [0, n)
            for (int i = 0; i + sl + 1 < tour.size(); i++) {
                int x1 = tour[i];
                int x2 = tour[i + 1];
                int y1 = tour[i + sl];
                int y2 = tour[i + sl + 1];
                for (int k = i + sl + 1; k + 1 < tour.size(); k++) {
                    if (k == i + sl + 1) {
                        ; //cout << "i: " <<  i << " k: " << k << endl; 
                    }
                    int z1 = tour[k];
                    int z2 = tour[k + 1];

                    int gain = dm[x1][x2] + dm[y1][y2] + dm[z1][z2];
                    int loss = dm[x1][y2] + dm[z1][x2] + dm[y1][z2];

                    int delta = loss - gain;
                    if (delta < 0) {

                        #if DEBUG == 1
                        double old_cost = tour_distance(tour, dm);
                        #endif

                        improved = true;
                        vector<int> saved_next_segment;
                        if (k + sl >= n) {break;}
                        for (int q = 0; q < sl; q++) {
                            //cout << "q: " << q << endl;
                            saved_next_segment.push_back(tour[q + k]);
                            tour[i + q + sl] = tour[i + q + 1];
                            tour[i + q + 1] = saved_next_segment[q];

                        }

                        saved_next_segment.clear();

                        #if DEBUG == 1
                        double new_cost = tour_distance(tour, dm);
                        cout << new_cost << " <-- " << old_cost << endl;
                        assert(new_cost < old_cost);
                        #endif

                        capital_delta += delta;
                    }
                }
            }
        }
    }
    return capital_delta;
}



int local_two_opt(vector<int>& tour, const vector<vector<int>>& dm, int parts) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;
    while (improved) {
        improved = false;
        for (int k = 0; k < parts + 1; k++) {
            for (int i = k*tour.size()/parts; i < (k+1)*tour.size()/parts - 1; i++) {
                for (int j = i + 1; j < (k+1)*tour.size()/parts - 1; j++) {

                    if (j >= tour.size() - 1 || i >= tour.size() - 1) {break;}
                    if (i <= 1 || j <= 1) {break;}
                    int i_m = mod((i - 1), n);
                    int j_p = mod((j + 1), n);  

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
                        //i -= 1;
                        //break;
                    }
                }
            }
        }
    }
    return capital_delta;
}

//Calculates change wrong
int three_opt(vector<int>& tour, const vector<vector<int>>& dm) {
    int n = tour.size();

    bool improved = true;
    int capital_delta = 0;

    while (improved) {
        improved = false;
        for (int i = 1; i < n ; ++i) {
            for (int j = i + 2; j < n - 2; ++j) {
                for (int k = j + 2; k < n - 2; ++k) {

                    int a, b, c, d, e, f;
                    a = tour[i - 1];
                    b = tour[i];
                    c = tour[j - 1];
                    d = tour[j];
                    e = tour[k - 1];
                    f = tour[k];
                    //cout << "a ... f" << a << " " << b << " " << c << " " << d << " " << e << " " << f  << endl;

                    int d0 = dm[a][b] + dm[c][d] + dm[e][f];
                    int d1 = dm[a][c] + dm[b][d] + dm[e][f];
                    int d2 = dm[a][b] + dm[c][e] + dm[d][f];
                    int d3 = dm[a][d] + dm[e][b] + dm[c][f];
                    int d4 = dm[f][b] + dm[c][d] + dm[e][a];

                    int delta = 0;

                    #if DEBUG == 1
                    double old_cost = tour_distance(tour, dm);
                    //cout << d0 << " " << d1 << " " << d2 << " " << d3 << " " << d4 << endl;
                    #endif

                    if (d0 > d1) {
                        std::reverse(tour.begin() + i, tour.begin() + j );
                        delta = -d0 + d1;
                    }
                    else if (d0 > d2) {
                        std::reverse(tour.begin() + j, tour.begin() + k );
                        delta = -d0 + d2;
                    }
                    
                    else if (d0 > d3) {
                        vector<int> t;
                        t.clear();
                        for (int iter = j; iter < k; iter++) {
                            t.push_back(tour[iter]);
                        }
                        for (int iter = i; iter < j; iter++) {
                            t.push_back(tour[iter]);
                        }

                        int iter = i;
                        for (auto& n: t) {
                            tour[iter] = n;
                            iter++;
                        }

                        delta = -d0 + d3;

                        #if DEBUG == 1
                        assert(all_different(tour));
                        #endif
                    }
                    
                    else if (d0 > d4) {
                        std::reverse(tour.begin() + i, tour.begin() + k );
                        delta = -d0 + d4;
                    }
                
                    if (delta < 0) {



                        //std::reverse(tour.begin() + i, tour.begin() + j + 1);
                        improved = true;

                        #if DEBUG == 1
                        double new_cost = tour_distance(tour, dm);
                        cout << old_cost << " -> " << new_cost << endl;
                        assert(new_cost < old_cost);
                        assert(new_cost == old_cost + delta);
                        #endif 
                        
                        capital_delta += delta;
    }   }   }   }   }

    return capital_delta;
}



int two_opt_move(vector<int>& tour, vector<vector<int>> dm, int i, int j) {


    #if DEBUG == 1
    vector<int> saved_tour(tour);
    #endif

    int n = tour.size();
    
    int i_m = mod((i - 1), n);
    int j_p = mod((j + 1), n);  

    int cost_1 = dm[tour[i_m]][tour[i]];
    int cost_2 = dm[tour[j]][tour[j_p]];

    int new_cost_1 = dm[tour[i]][tour[j_p]];
    int new_cost_2 = dm[tour[i_m]][tour[j]];

    int delta = (new_cost_1 + new_cost_2) - (cost_1 + cost_2);

    std::reverse(tour.begin() + i , tour.begin() + j + 1);

    #if DEBUG == 1
    assert(tour_distance(tour, dm) ==  delta + tour_distance(saved_tour, dm));
    if (delta < 0) assert(tour_distance(tour, dm) < tour_distance(saved_tour, dm));
    #endif

    return delta;
}



int two_opt_move(tour_o& tour, vector<vector<int>> dm, int i, int j) {


    #if DEBUG == 1
    tour_o saved_tour(tour);
    #endif

    int n = tour.size();
    
    int i_m = mod((i - 1), n);
    int j_p = mod((j + 1), n);  

    #if DEBUG == 1
    assert(n == dm.size());
    assert(i_m >= 0 && i_m < n);
    assert(j_p < n && j_p >= 0);
    assert(tour.at(i_m) >= 0);
    assert(tour.at(j_p) < n);
    /*
    cout << i_m << endl;
    cout << tour.at(i_m) << endl;
    cout << j_p << endl;
    cout << tour.at(j_p) << endl;
    cout << i << endl;
    cout << tour.at(i) << endl;
    cout << j << endl;
    cout << tour.at(j) << endl;
    */
    #endif

    int cost_1 = dm[tour.at(i_m)][tour.at(i)];
    int cost_2 = dm[tour.at(j)][tour.at(j_p)];

    int new_cost_1 = dm[tour.at(i)][tour.at(j_p)];
    int new_cost_2 = dm[tour.at(i_m)][tour.at(j)];

    int delta = (new_cost_1 + new_cost_2) - (cost_1 + cost_2);

    tour.reverse(i, j);

    #if DEBUG == 1
    assert(tour_distance(tour, dm) ==  delta + tour_distance(saved_tour, dm));
    if (delta < 0) assert(tour_distance(tour, dm) < tour_distance(saved_tour, dm));
    #endif

    return delta;
}

vector<int> k_opt(vector<int> tour, const vector<vector<int>>& dm, const vector<vector<int>>& nm, const vector<set<int>>& ns, int depth) {
    // 4 3 2 scheme

    int ts = tour.size();
    assert(ts >= 1);
    vector<int> new_tour(tour);
    int delta = 0;
    int a;
    int b;

    for (int i = 0; i < depth; i++) {
        //a = b;
        b = ts - 1;
        //TOFIX: allow a == 0

        a = 1 + (int) ((ts - 1) * (rand() / (RAND_MAX + 1.0)));
        //a = 0 + (int) ((nbors - 1) * (rand() / (RAND_MAX + 1.0)));
        //a = nm[b][a] + 1;
        //a = rand() %
        delta += two_opt_move(new_tour, dm, a, b);
        if (delta < 0) {return tour;}
    }

    #if VERBOSE == 1
    auto t1 = chrono::high_resolution_clock::now();
    #endif

    //delta += two_opt_nn(new_tour, dm, nm, ns);
    delta += two_opt(new_tour, dm);
    delta += node_shift(new_tour, dm);
    /*
    if (tour.size() < 200) {
        delta += three_opt(new_tour, dm);
    }
    */

    #if VERBOSE == 1
    auto t2 = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time_since = t2 - t1;
    cout << "two opt took " << time_since.count() << endl;
    #endif

    //cout << tour_distance(tour, dm) << " " <<  tour_distance(saved_tour, dm) << " " << delta;
    //cout << "Kopt delta " << delta << endl;

    if (delta < 0) {
        #if DEBUG == 1
        vector<int> saved_tour(tour);
        #endif
        tour = new_tour;
        #if DEBUG == 1
        assert(tour_distance(tour, dm) < tour_distance(saved_tour, dm));
        #endif
    }
    
    return tour;
}


tour_o k_opt(tour_o tour, const vector<vector<int>>& dm, const vector<vector<int>>& nm, const vector<set<int>>& ns, int depth) {
    // 4 3 2 scheme

    int ts = tour.size();
    assert(ts >= 1);
    tour_o new_tour(tour.id);
    int delta = 0;
    int a;
    int b;

    for (int i = 0; i < depth; i++) {
        //a = b;
        b = ts - 1;
        //TOFIX: allow a == 0
        a = 1 + (int) ((ts - 1) * (rand() / (RAND_MAX + 1.0)));
        //a = rand() %
        delta += two_opt_move(new_tour, dm, a, b);
        if (delta < 0) {return tour;}
    }

    #if VERBOSE == 1
    auto t1 = chrono::high_resolution_clock::now();
    #endif

    //cout << tour_distance(tour, dm) << " " <<  tour_distance(new_tour, dm) << " " << delta << endl;
    //cout << "Kopt delta " << delta << endl;

    //delta += two_opt_nn(new_tour, dm, nm, ns);
    delta += two_opt(new_tour, dm);

    #if VERBOSE == 1
    auto t2 = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> time_since = t2 - t1;
    cout << "two opt took " << time_since.count() << endl;
    #endif

    //cout << tour_distance(tour, dm) << " " <<  tour_distance(new_tour, dm) << " " << delta << endl;
    //cout << "Kopt delta " << delta << endl;

    if (delta < 0) {

        #if DEBUG == 1
        int old_cost = tour_distance(tour, dm);
        #endif
        
        tour = new_tour;

        #if DEBUG == 1
        assert(tour_distance(new_tour, dm) < old_cost);
        #endif
    }
    
    return tour;
}



vector<int> anneal(vector<int> tour, vector<vector<int>> dm, double time_limit) {

    #if VERBOSE == 1
    int counter = 0;
    #endif

    auto start_time = chrono::high_resolution_clock::now();

    int n = tour.size();
    int cost = tour_distance(tour, dm);

    double cooling_rate = 0.4;
    double temperature  = 200;

    chrono::duration<double, milli> time_since;

    while (true) {

        auto curr = chrono::high_resolution_clock::now();
        time_since = curr - start_time;
        if ((time_since).count() > time_limit) {break;}

        int v1 = 1 + rand() % (n - 1);
        int u1 = 1 + rand() % (n - 1);
        v1 = min(v1, u1);
        u1 = max(v1, u1);

        //Random swap of 2 elements
        int gain_rand = dm[tour[v1 - 1]][tour[v1]] + dm[tour[v1]][tour[v1 + 1]];
        gain_rand    += dm[tour[u1 - 1]][tour[u1]] + dm[tour[u1]][tour[u1 + 1]];

        int loss_rand = dm[tour[v1]][tour[u1 - 1]] + dm[tour[v1]][tour[u1 + 1]];
        loss_rand    += dm[tour[u1]][tour[v1 - 1]] + dm[tour[u1]][tour[v1 + 1]];

        int v2 = 1 + rand() % (n - 2);
        int u2 = 1 + rand() % (n - 2);

        v2 = min(v2, u2);
        u2 = max(v2, u2);

        //TOFIX
        if (u2 == v2 || u1 == v1) {continue;}

        int gain_2opt  = dm[tour[v2 - 1]][tour[v2]];
        gain_2opt     += dm[tour[u2]][tour[u2 + 1]];

        int loss_2opt  = dm[tour[v2]][tour[u2 + 1]];
        loss_2opt     += dm[tour[v2 - 1]][tour[u2]];

        int new_rand = cost + gain_rand - loss_rand;
        double accept_p_rand = exp(-(new_rand - cost) / temperature);

        /*
        if (accept_p_rand > rand()/RAND_MAX || new_rand < cost) {
            //update
            cost = new_rand;
            int s = tour[u1];
            tour[v1] = tour[u1];
            tour[u1] = s;
            cost = cost - gain_rand + loss_rand;

            //assert(tour_distance(tour, dm) == cost);
            continue;
        }            cost = tour_distance(tour, dm);
        */
        

        int new_2opt = cost - gain_2opt + loss_2opt;
        double accept_p_2opt = exp(-(new_2opt - cost) / temperature);

        if (accept_p_2opt > rand()/RAND_MAX || new_2opt < cost) {
            //update
            reverse(tour.begin() + v2, tour.begin() + u2 + 1);
            cost = new_2opt;
        }   

        temperature = max(0., temperature - cooling_rate);

        #if DEBUG == 1
        assert(tour_distance(tour, dm) == cost);
        #endif

        #if VERBOSE == 1
        if (counter % 1000 == 0) {
            cout << "cost at time c: " << counter << " is: " << cost << endl;
            cout << "temperature:" << temperature << endl;
            cout << "time: " << time_since.count() << endl;
        }
        counter++;
        #endif
    }

    #if VERBOSE == 1
    cout << "anneal ran n times: " << counter << endl;
    cout << "temperature at end" << temperature << endl;
    cout << "in time " << time_since.count() << endl;
    cout << "new cost " << tour_distance(tour, dm) << endl;
    #endif

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



    static vector<vector<int>> dm = create_distance_matrix(points);
    static vector<vector<int>> nm = create_nearest_neighbor_matrix(dm, 10);
 
    static vector<int> tour(n);

    for (int i = 0; i < n; i++) {
        tour[i] = i;
    }

    //TOFIX constructor

    tour = greedy_tour(dm);

    two_opt(tour, dm);
    /*
    if (n < 200) {
        three_opt(tour, dm);
    }
    */
    node_shift(tour, dm);



    //greedy_tour(&tour, dm);

    #if DEBUG == 1
    assert(tour_.id.size() == n);
    assert(tour_.pos.size() == n);
    #endif

    auto t1 = chrono::high_resolution_clock::now();
    auto t2 = chrono::high_resolution_clock::now();

    #if VERBOSE == 1
    int best_tour_dist = tour_distance(tour, dm);//*max_element(dm[0].begin(), dm[0].end()) * n  + 1;
    #endif


    vector<set<int>> ns;
    for (int i = 0; i < n; i++){
        ns.push_back(set<int>());
        for (auto& n: nm[i]) {
            ns[i].insert(n);
        }
    }

    int depth = 3;
    int i = 0;

    #if VERBOSE == 1
    cout << "Start optimize, cost is " << tour_distance(tour, dm) << endl;
    cout << "(Cost of normal tour was " << tour_distance(tour_, dm) << endl;
    #endif

    while (true) {

        tour = k_opt(tour, dm, nm, ns, depth);
        //tour = anneal(tour, dm, 1700);
        //two_opt(tour, dm);

        t2 = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> time_since = t2 - t1;

        #if VERBOSE == 1
        int new_dist = tour_distance(tour, dm);
        if (new_dist < best_tour_dist) {
            cout << new_dist << " <- " << best_tour_dist << " in " << time_since.count() << endl;
            best_tour_dist = new_dist;
        }
        #endif

        //or_opt(tour, dm);
        if (time_since.count() > 1890) {
            break;
        }
        #if VERBOSE == 1
        i++;
        #endif
    }

    two_opt(tour, dm);

    #if VERBOSE == 1
    cout << tour_distance(tour, dm) << endl;
    cout << "n steps " << i << endl;
    #endif

    for (int i = 0; i < n; i++) {
        cout << tour.at(i);
        if (i < n-1)
            cout << "\n";
    }

  return 0;

}
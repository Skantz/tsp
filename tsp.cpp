#include <iostream>
#include <sstream>
#include <istream>
#include <vector>
#include <math.h> 
#include <algorithm> 
#include <assert.h>

using namespace std;

#define DEBUG 1

//x1, y1, x2, y2..
vector<vector<int>> create_distance_matrix(vector<double> points) {
    int n = points.size() / 2;
    vector<vector<int>> dm;
    for (int i = 0; i < points.size(); i++) {
        dm.push_back(vector<int>());
    }
    for (int i = 0; i < points.size(); i++) {
        for (int j = 0; j < points.size(); j++) {
            double x_dist = (points[i*2] - points[j*2]);
            double y_dist = (points[i*2 + 1] - points[j*2 + 1]);
            dm[i][j] = pow(pow(x_dist, 2) + pow(y_dist, 2), 0.5);
        }
    }
    return dm;
}

int tour_distance(vector<int> tour, vector<vector<int>> dm) {
    
    int s = 0;
    for (int i = 0; i < tour.size() - 1; i++) {
        s += dm[tour[i]][tour[i + 1]];
    }
    return s + dm[tour[tour.size() - 1]][tour[0]];
}

void two_opt(vector<int>& tour, vector<vector<int>>& dm) {
    bool improved;
    for (int i = 0; i < tour.size(); i++) {
        for (int j = 0; j < tour.size(); j++) {
            int cost_1 = dm[tour[i - 1]][tour[i]];
            int cost_2 = dm[tour[j]][tour[j + 1]];

            int new_cost_1 = dm[tour[i]][tour[j + 1]];
            int new_cost_2 = dm[tour[i-1]][tour[j]];

            if ((new_cost_1 + new_cost_2) < (cost_1 + cost_2)) {
                //SWAP
                #if DEBUG == 1
                double old_cost = tour_distance(tour, dm);
                #endif
                std::reverse(tour.begin() + i +1, tour.begin() + j);
                improved = true;
                #if DEBUG == 1
                double new_cost = tour_distance(tour, dm);
                assert(new_cost < old_cost);
                #endif 
            }
        }
    }
    return;
}


int main(int argc, char *argv[]) {

    vector<double> points;

    int n;
    cin >> n;

    double val;
    for (int i = 0; i < 2*n; i++) {
        std::cin >> val;
        points.push_back(val);
    }

    vector<vector<int>> dm = create_distance_matrix(points);
    
    vector<int> tour;
    for (int i = 0; i < n; i++) {
        tour.push_back(i);
    }


    if (argc > 1) {
        cout << tour_distance(tour, dm);
    }

    else {
        for (int i = 0; i < n; i++) {
            cout << tour[i];
            if (i < n-1)
                cout << "\n";
            }
    }

  return 0;

}

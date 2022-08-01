#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
using namespace std;
using namespace std::chrono;
int seed = 1322;
// random_device rd;
// mt19937 generator(rd());
default_random_engine generator (seed);
uniform_real_distribution<float> distribution(0.0,1.0);

class DisjointSet {
  public:
    int *rank, *parent, N;
    DisjointSet(int n){
      rank = new int[n];
      parent = new int[n];
      N = n;
      makeSet();
    }

    //MAKESET function
    void makeSet(){
      for (int i=0; i < N; i++){
        parent[i] = i;
        rank[i] = 0;
      }
    }

    //FIND function
    int find(int k){
      if (parent[k] != k){
        parent[k] = find(parent[k]);
      }
      return parent[k];
    }

    //UNION function
    void Union(int x, int y){
      int xp = find(x);
      int yp = find(y);

      if (xp == yp){
        return;
      }
      if (rank[xp] < rank[yp]){
        parent[xp] = yp;
      } else if (rank[xp] > rank[yp]){
        parent[yp] = xp;
      } else {
        parent[yp] = xp;
        rank[xp] += 1;
      }
    }
};

class Graph {
  public:
    int V; //number of vertices
    vector<vector<float> > edges; //vector of edges
    float maxw=0; //max weighted edge used in MST
    Graph(int vertices) {
      V = vertices;
    }

    //add edge to the graph -- u and v are 2 endpoints, w is weight
    void addEdge(float u, float v, float w){
      edges.push_back({w,u,v});
    }

    float kruskal(){
      float cost = 0;

      DisjointSet S(V);

      sort(edges.begin(), edges.end());

      for (auto edge : edges){
        float w = edge[0];
        float u = edge[1];
        float v = edge[2];

        if (S.find(u) != S.find(v)){
          S.Union(u,v);
          cost += w;
          if (w > maxw){
            maxw = w;
          }
        }
      }
      return cost;
    }
};

//generate 0d graph nodes and weights
void gen0d(Graph& G, int n){
  for (int i=0; i < n-1; i++){
    for (int j=i+1; j < n; j++){
      float weight = distribution(generator);
      G.addEdge(i,j,weight);
      cout << i << "-" << j << " " << weight << endl; 
    }
  }
}

//generate 2d graph nodes and weights
void gen2d(Graph& G, int n){
  vector<vector<float> > points; //vector of points
  for (int i=0; i < n; i++){
    float x = distribution(generator);
    float y = distribution(generator);
    points.push_back({x,y});
  }
  for (int i=0; i < n-1; i++){
    for (int j=i+1; j < n; j++){
      float weight = sqrt(pow((points[i][0]-points[j][0]),2) + pow((points[i][1]-points[j][1]),2));
      G.addEdge(i,j,weight);
    }
  }
}

//generate 3d graph nodes and weights
void gen3d(Graph& G, int n){
  vector<vector<float> > points; //vector of points
  for (int i=0; i < n; i++){
    float x = distribution(generator);
    float y = distribution(generator);
    float z = distribution(generator);
    points.push_back({x,y,z});
  }
  for (int i=0; i < n-1; i++){
    for (int j=i+1; j < n; j++){
      float weight = sqrt(pow((points[i][0]-points[j][0]),2) + pow((points[i][1]-points[j][1]),2) + pow((points[i][2]-points[j][2]),2));
      G.addEdge(i,j,weight);
    }
  }
}

//generate 4d graph nodes and weights
void gen4d(Graph& G, int n){
  vector<vector<float> > points; //vector of points
  for (int i=0; i < n; i++){
    float x = distribution(generator);
    float y = distribution(generator);
    float z = distribution(generator);
    float t = distribution(generator);
    points.push_back({x,y,z,t});
  }
  for (int i=0; i < n-1; i++){
    for (int j=i+1; j < n; j++){
      float weight = sqrt(pow((points[i][0]-points[j][0]),2) + pow((points[i][1]-points[j][1]),2) + pow((points[i][2]-points[j][2]),2) + pow((points[i][3]-points[j][3]),2));
      G.addEdge(i,j,weight);
    }
  }
}

int main(int argc, char* argv[]){
  //check arguments
  if (argc != 5){
    cout << "Invalid arguments, should be './randmst 0 numpoints numtrials dimension'" << endl;
    return 0;
  }

  int flag = atoi(argv[1]);
  int numpoints = atoi(argv[2]);
  int numtrials = atoi(argv[3]);
  int dimension = atoi(argv[4]);

  //total cost of MST through all trials
  float mst = 0;

  //start timer
  auto start = high_resolution_clock::now();

  //interate through trials
  for (int t=0; t < numtrials; t++){
    //instantiate a graph
    Graph G(numpoints);

    //create nodes and edges based on dimension
    if (dimension == 0){
      gen0d(G,numpoints);
    } else if (dimension == 2){
      gen2d(G,numpoints);
    } else if (dimension == 3){
      gen3d(G,numpoints);
    } else {
      gen4d(G,numpoints);
    }

    mst += G.kruskal();
    cout << "max edge weight = " << G.maxw << endl;
  }

  //stop timer
  auto stop = high_resolution_clock::now();

  //calculate time executed
  auto duration = duration_cast<milliseconds>(stop - start);

  //average MST cost
  float avgMST = mst/numtrials;

  cout << avgMST << " " << numpoints << " " << numtrials << " " << dimension << endl;
  cout << "Time taken: "<< duration.count() << " milliseconds" << endl;
  return 0;
}

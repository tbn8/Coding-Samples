#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include <chrono>
using namespace std;
using namespace std::chrono;
float inf = std::numeric_limits<float>::max();

//*******************************
//Uniform real random generator
//*******************************
//set seed for code testing only
// int seed = 1322;
// default_random_engine generator(seed);

//generate random seed number for every trial
random_device rd;
mt19937 generator(rd());

//generate uniform real distribution
uniform_real_distribution<float> distribution(0.0,1.0);


//*******************************
//Graph structures and functions
//*******************************
//node in adjacency list
struct Node {
  int dest; //destination node
  float weight; //weight of edge from source node to destination node
  Node *next; //pointer to next node
};

struct AdjacencyList {
  Node *head; //pointer to head of list
};

//graph is an array of adjacency list
struct Graph {
  int V; //number of vertices
  AdjacencyList *array;
};

//add a node to adjacency list
Node *addNode(int dest, float weight){
  Node *node = (Node*) malloc(sizeof(Node));
  node->dest = dest;
  node->weight = weight;
  node->next = NULL;
  return node;
}

//build graph
Graph *buildGraph(int V){
  Graph *G = (Graph*) malloc(sizeof(Graph));
  G->V = V;
  G->array = (AdjacencyList*) malloc(V * sizeof(AdjacencyList));
  for (int i=0; i < V; i++){
    G->array[i].head = NULL;
  }
  return G;
}

//add an edge from u to v with weight w to graph
void addEdge(Graph *G, int u, int v, float w){
  //add node v to the top of adjacency list of u
  Node *node = addNode(v,w);
  node->next = G->array[u].head;
  G->array[u].head = node;

  //add node u to the top of adjacency list of v
  node = addNode(u,w);
  node->next = G->array[v].head;
  G->array[v].head = node;
}

//*****************************************
//Min binary heap structures and functions
//*****************************************
struct HeapNode {
  int vertex; //represents vertex in graph
  float key; //represents min distance from root of MST
};

struct MinHeap {
  int capacity; //max size of heap
  int size; //current number of nodes in heap
  int *loc; //pointer to node location
  HeapNode **heapArr; //pointer to heap array
};

//insert a node to the heap
HeapNode *insertHeap(int v, float key){
  HeapNode *node = (HeapNode*) malloc(sizeof(HeapNode));
  node->vertex = v;
  node->key = key;
  return node;
}

//build min heap
MinHeap *buildHeap(int capacity){
  MinHeap *heap = (MinHeap*) malloc(sizeof(MinHeap));
  heap->capacity = capacity;
  heap->size = 0;
  heap->loc = (int*) malloc(capacity * sizeof(int));
  heap->heapArr = (HeapNode**) malloc(capacity * sizeof(HeapNode*));
  return heap;
}

//rearrange nodes in heap
void minHeapify(MinHeap *heap, int i){
  int min = i; //set min node to current index
  int left = i*2 + 1; //left child of current index
  int right = i*2 + 2; //right child of current index

  //compare index to left or right child
  if (left < heap->size && heap->heapArr[left]->key < heap->heapArr[min]->key){
    min = left;
  }
  if (right < heap->size && heap->heapArr[right]->key < heap->heapArr[min]->key){
    min = right;
  }

  //swap min nodes
  if (min != i) {
    HeapNode *minNode = heap->heapArr[min];
    HeapNode *indexNode = heap->heapArr[i];
    heap->loc[minNode->vertex] = i;
    heap->loc[indexNode->vertex] = min;
    heap->heapArr[min] = indexNode;
    heap->heapArr[i] = minNode;
    minHeapify(heap, min);
  }
}

//extract min node from heap
HeapNode *extractMin(MinHeap *heap){
  if (heap->size == 0){
    return NULL;
  }

  //store min node
  HeapNode *root = heap->heapArr[0];

  //swap root and last
  HeapNode *last = heap->heapArr[heap->size - 1];
  heap->heapArr[0] = last;
  heap->loc[root->vertex] = heap->size - 1;
  heap->loc[last->vertex] = 0;

  //reduce heap size and heapify
  heap->size --;
  minHeapify(heap,0);

  return root;
}

//bubble up a node
void bubbleUp(MinHeap *heap, int v, float key){
  int i = heap->loc[v]; //index of vertex in heap array
  int parent = (i-1)/2; //index of parent
  heap->heapArr[i]->key = key;

  //swap node upwards
  while (i>0 && heap->heapArr[i]->key < heap->heapArr[parent]->key){
    heap->loc[heap->heapArr[i]->vertex] = parent;
    heap->loc[heap->heapArr[parent]->vertex] = i;
    HeapNode *temp = heap->heapArr[i];
    heap->heapArr[i] = heap->heapArr[parent];
    heap->heapArr[parent] = temp;
    i = parent;
    parent = (i-1)/2;
  }
}

//*****************
//Prim's Algorithm
//*****************
float Prim(Graph *G){
  int V = G->V; //number of vertices
  int parent[V]; //array of MST
  float key[V]; //key values to store distance from root MST
  MinHeap *heap = buildHeap(V);
  heap->size = V;

  //add first node to heap, setting key aka distance to itself to 0
  key[0] = 0;
  heap->heapArr[0] = insertHeap(0, key[0]);
  heap->loc[0] = 0;

  //add all other nodes to heap, setting key to infinity
  for (int v=1; v < V; v++){
    parent[v] = -1;
    key[v] = inf;
    heap->heapArr[v] = insertHeap(v, key[v]);
    heap->loc[v] = v;
  }

  while (heap->size != 0){
    HeapNode *minNode = extractMin(heap);
    int u = minNode->vertex;

    //go through adjacent nodes of u and update keys
    Node *node = G->array[u].head;
    while (node != NULL) {
      int v = node->dest;

      //if v is still in heap and weight of u->v < key of v then update
      if (heap->loc[v] < heap->size && node->weight < key[v]){
        key[v] = node->weight;
        parent[v] = u;
        bubbleUp(heap, v, key[v]);
      }

      node = node->next;
    }
  }

  //add up MST cost
  float cost = 0;
  // float maxw = 0;
  for (int i=1; i < V; i++){
    cost += key[i];
    // if (key[i] > maxw){
    //   maxw = key[i];
    // }
  }
  // cout << "max weight used = " << maxw << endl;
  return cost;

}

//****************************************
//Generate random graph edges and weights
//****************************************
//generate cutoff value to prune edges
float genK(int n, int dim, int flag){
  if (flag == 0 || n < 128){
    return 1.0;
  }
  if (dim == 0){
    return 1/pow(n,0.59);
  }
  if (dim == 2){
    return 1/pow(log2(pow(n,0.25)),2.5);
  }
  if (dim == 3){
    return 1/(log2(n)-1.6*log2(log2(n)));
  }
  if (dim == 4){
    return 1/log2(sqrt(n/12));
  }
  return 1.0;
}

//generate 0d graph edges and weights
void gen0d(Graph *G, int n, float k){
  for (int i=0; i < n-1; i++){
    for (int j=i+1; j < n; j++){
      float weight = distribution(generator);
      if (weight < k){
        addEdge(G,i,j,weight);
      }
    }
  }
}

//generate 2d graph edges and weights
void gen2d(Graph *G, int n, float k){
  vector<vector<float> > points; //vector of points
  for (int i=0; i < n; i++){
    float x = distribution(generator);
    float y = distribution(generator);
    points.push_back({x,y});
  }
  for (int i=0; i < n-1; i++){
    for (int j=i+1; j < n; j++){
      float weight = sqrt(pow((points[i][0]-points[j][0]),2) + pow((points[i][1]-points[j][1]),2));
      if (weight < k){
        addEdge(G,i,j,weight);
      }
    }
  }
}

//generate 3d graph edges and weights
void gen3d(Graph *G, int n, float k){
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
      if (weight < k){
        addEdge(G,i,j,weight);
      }
    }
  }
}

//generate 4d graph edges and weights
void gen4d(Graph *G, int n, float k){
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
      if (weight < k){
        addEdge(G,i,j,weight);
      }
    }
  }
}

//***************
//MAIN EXECUTION
//***************
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

  //cutoff value for edge pruning
  float k = genK(numpoints, dimension, flag);
  cout << "k = " << k << endl;

  //interate through trials
  for (int t=0; t < numtrials; t++){
    //instantiate a graph
    Graph *G = buildGraph(numpoints);

    //create nodes and edges based on dimension
    if (dimension == 0){
      gen0d(G, numpoints, k);
    } else if (dimension == 2){
      gen2d(G, numpoints, k);
    } else if (dimension == 3){
      gen3d(G, numpoints, k);
    } else {
      gen4d(G, numpoints, k);
    }

    mst += Prim(G);
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

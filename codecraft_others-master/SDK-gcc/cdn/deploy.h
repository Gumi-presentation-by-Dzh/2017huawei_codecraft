#ifndef __ROUTE_H__
#define __ROUTE_H__

#include "lib_io.h"
#include <vector>

using namespace std;

const int MAX_V = 2010;
const int inf = 0x3f3f3f3f;  // I think this variable means infinite
const long long infll = 0x3f3f3f3f3f3f3f3f;
const double c1 = 1.0;
const double c2 = 1.6;
const double w = 0.9;

/* @v is neighbor's id, @c is link's cost */
struct EdgeInfo {
    int v;
    int c;
    EdgeInfo(int _v, int _c): v(_v), c(_c) {}
};

// @v is the corresponding network node, @w is the need
struct CustomerNodeInfo {
    int v;
    int w;
    CustomerNodeInfo(int _v, int _w): v(_v), w(_w) {}
};

struct Edge {
    int t;          // neighbor id
    int u;          // bandwidth
    int c;          // cost per Gbps 
    int U;          // bandwidth
    int C;          // cost per Gbps
    Edge *next;     // the next edge of this vertex
    Edge *pair;     // the reverse edge of the edge, the cost is negative-@c, the bandwidth sum is zero
};

class Fuck {
public:
    void add_server(vector<int> & v);
    void add_edge(int u, int v, int w, int c);
    long long costflow();
    void print_flow(vector<vector<int> > &node, vector<int> &flow);
    void readtopo(char * topo[MAX_EDGE_NUM], int line_num);
    void spfa();
    void kmeans(int k, vector<int> & clusters);

    int need_flow;       // total need of video resources
    int node_num;        // network node num
    int edge_num;        // network edge num
    int customer_num;    // customer node num
    long long server_cost;  // constant cost of building server
private:
    int aug(int u, int m);
    bool modlabel();

    vector<vector<EdgeInfo> > graph;

    vector<CustomerNodeInfo> customer_nodes; // all customer nodes' infomations

    vector<vector<int> > d;
    int vis[MAX_V]; // used in spfa(), see the usage in spfa()'s definition

    // store the real edge info in @epool, @e only contains the pointer to the correct edges
    // the indices bigger than node_num 
    Edge epool[MAX_EDGE_NUM * 5];
    Edge *e[MAX_V]; // every vertex's induced edge
    
    int psz; // the current index in @epool
    int psz_tmp;  // store the @psz temporary

    // @s is the index of super source, and @t is super sink
    int s;   // total node num, including network nodes and customer nodes
    int t;   // t's init value = s+1 == node_num + customer_num + 1
    
    int flow;
    int dist;
    int D[MAX_V];
    long long cost;
};

class Particle {
public:
    Particle(int length=0);
    Particle(int length, vector<int> & vi, Fuck* & fuck);

    vector<double> v;       // 
    vector<double> v_best;  // this particle's local best position dimensions
    vector<double> vp;      // I think this is a temproary position

    long long cost_best;
    long long cost;
};

class XJBS {
public:
    XJBS(Fuck & fuck);
    void get_best(vector<int> & server);
    void addone(vector<int> & v);
    void run1();
    void run2();
    double initial(int size);
    inline void reproduction();
    vector<Particle> p;

private:
    void decode(vector<double> & vd, vector<int> & vi);

    void GA_cross(Particle & s1, Particle & s2);
    void OBMA(Particle & s);
    inline void PSO_update(Particle & s);
    inline void updateone(Particle & s);


    Particle gbest;     // global best position
    int l;              // the @fuck's nodenum
    int max_p_size;     // equals 10 during running and not be changed
    int cnt;
    double PSO_c1;      // PSO formula's first accelerate rate
    double PSO_c2;      // PSO formula's second accelerate rate
    double PSO_w;       // PSO formula's weight of current speed inertia
    Fuck *fuck;
};

template <class T>
void knuth_shuffle(vector<T> & v);

bool cmp(const Particle & p1, const Particle & p2);

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename);

#endif


#ifndef __ROUTE_H__
#define __ROUTE_H__

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <string>
#include <queue>
#include <deque>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <climits>
#include "lib_io.h"
#include "lib_time.h"
#include "deploy.h"

#define INF 0x3f3f3f3f
#define MAX_ITERATION 100000

using namespace std;

const int MAX_V = 2010;

struct edge;
struct consumer;
struct graph;
struct costflow;
struct PSO;
struct Particle;

typedef int vid;
typedef long long LL;

void deploy_server(char * graph[MAX_EDGE_NUM], int edge_num, char * filename);
/* these functions are used for debugging */
template <typename T>
void pvector(vector<T> &v, int fd = 1);


struct edge
{
    edge(vid _v, int _cap, int _c, int _e = 0, int _bw = 0, int _cost = INF):
        v(_v), cap(_cap), c(_c), re(_e), bw(_bw), cost(_cost){
            ifre = false;
        }

    vid v;
    int cap;    // the current capacity(used in zkw)
    int c;      // the current cost with update(used in zkw)
    int bw;     // total bandwidth
    int cost;   // cost per Gbps
    int re;   // the pair edge with same ends
    bool ifre;  // indicate whether this is a reverse edge
};

struct consumer
{
    consumer():v(-1), need(INF) {}
    consumer(vid _v, int _need):v(_v), need(_need){}
    vid v;
    int need;
};

struct graph
{
    void init(char *topo[MAX_EDGE_NUM], int linenum);
    void addEdge(vid u, vid v, int bw, int cost, int flag);
    void cluster(int k, vector<vid> &cores);
    void allpairs_sp();
    void make_supersource(vector<vid> &sours);
    void make_supersink();

    // @n_node_num front nodes are network nodes, others are consumer nodes 
    vector<vector<edge>> e;
    vector<consumer> consumers;  // all consumer node infomations
    vector<vector<int>> d;

    int n_node_num;          // network node number
    int c_node_num;          // consumer node number
    int edge_num;
    int totalneed;
    int server_cost;         // constant cost of server

    vid ssour;               // super source
    vid sdest;               // super sink
    int have_superends;      // a flag
};

/* use "zkw" algorithm of Min-cost Max-flow problem */
struct costflow
{
    costflow(graph &_g);
    LL compute();
    LL getflowcost();
    void print_flow(vector<vector<vid>> &nodes, vector<int> &flow);
    bool modlabel();
    int augmenting(int u, int m);
    void resetvis(int value);
    void resetD(int value);

    graph &g;
    const vector<vid> servers;
    vector<int> vis;
    //vector<int> D;
    vector<int> D;
    vid dest;
    int res;
    int flow;
    int dist;
    private:
        LL cost;
};

struct Particle {
    Particle(int length=0);
    Particle(int length, vector<int> & vi, graph& g);

    vector<double> v;
    vector<double> v_best;
    vector<double> vp;

    long long cost_best;
    long long cost;
};

struct PSO
{
    PSO(graph &g_in);
    double init(int size); // initial(int size)
    void add_candidate(vector<int> &cores); // addone
    void phase(int flat);
    inline void reproduction();

    void get_best(vector<int> &servers);

    void decode(vector<double> & vd, vector<int> & vi);
    void GA_cross(Particle & s1, Particle & s2);
    void OBMA(Particle & s);
    inline void PSO_update(Particle & s);
    inline void updateone(Particle & s);


    graph &g;
    vector<Particle> candidates;///instead p;
    Particle gbest; 
    double c1, c2, w;   //pso
    double last_second;
    double first_second;
    // temp 
    int max_size;  // max_p_size
    int cnt;
    int l;
};

bool cmp(const Particle& p1, const Particle& p2);


#endif

